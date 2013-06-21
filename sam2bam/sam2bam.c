/* The MIT License

   Copyright (c) 2013 by Margus Lukk <margus.lukk@cruk.cam.ac.uk>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <pthread.h>
#include <assert.h>
#include "sam.h"
#include "bgzf.h"
#include "khash.h"
#include "bam_endian.h"
#include "kstring.h"

KHASH_SET_INIT_STR(rg)

// When counting records instead of printing them,
// data passed to the bam_fetch callback is encapsulated in this struct.

typedef struct {
  bam_header_t *header;
  int *count;
} count_func_data_t;

typedef khash_t(rg) *rghash_t;

// FIXME: we'd better use no global variables...

static rghash_t g_rghash = 0;
static int g_min_mapQ = 0, g_flag_on = 0, g_flag_off = 0;
static float g_subsam = -1;
static char *g_library, *g_rg;
static void *g_bed;

// following set of global variables is from bgzf.h
typedef int8_t bgzf_byte_t;

static const int DEFAULT_BLOCK_SIZE = 64 * 1024;
static const int MAX_BLOCK_SIZE = 64 * 1024;

static const int BLOCK_HEADER_LENGTH = 18;
static const int BLOCK_FOOTER_LENGTH = 8;

static const int GZIP_ID1 = 31;
static const int GZIP_ID2 = 139;
static const int CM_DEFLATE = 8;
static const int FLG_FEXTRA = 4;
static const int OS_UNKNOWN = 255;
static const int BGZF_ID1 = 66; // 'B'

static const int BGZF_ID2 = 67; // 'C'

static const int BGZF_LEN = 2;
static const int BGZF_XLEN = 6; // BGZF_LEN+4

static const int GZIP_WINDOW_BITS = -15; // no zlib header
static const int Z_DEFAULT_MEM_LEVEL = 8;



// following function is copied one to one from samtools
static inline int __g_skip_aln(const bam_header_t *h, const bam1_t *b)
{
  if (b->core.qual < g_min_mapQ || ((b->core.flag & g_flag_on) != g_flag_on) || (b->core.flag & g_flag_off))
    return 1;
  if (g_bed && b->core.tid >= 0 && !bed_overlap(g_bed, h->target_name[b->core.tid], b->core.pos, bam_calend(&b->core, bam1_cigar(b))))
    return 1;
  if (g_subsam > 0.) {
    int x = (int)(g_subsam + .499);
    uint32_t k = __ac_X31_hash_string(bam1_qname(b)) + x;
    if (k%1024 / 1024.0 >= g_subsam - x) return 1;
  }
  if (g_rg || g_rghash) {
    uint8_t *s = bam_aux_get(b, "RG");
    if (s) {
      if (g_rg) return (strcmp(g_rg, (char*)(s + 1)) == 0)? 0 : 1;
      if (g_rghash) {
	khint_t k = kh_get(rg, g_rghash, (char*)(s + 1));
	return (k != kh_end(g_rghash))? 0 : 1;
      }
    }
  }
  if (g_library) {
    const char *p = bam_get_library((bam_header_t*)h, b);
    return (p && strcmp(p, g_library) == 0)? 0 : 1;
  }
  return 0;
}


// Following function copied one to one from samtools
static char *drop_rg(char *hdtxt, rghash_t h, int *len)
{
  char *p = hdtxt, *q, *r, *s;
  kstring_t str;
  memset(&str, 0, sizeof(kstring_t));
  while (1) {
    int toprint = 0;
    q = strchr(p, '\n');
    if (q == 0) q = p + strlen(p);
    if (q - p < 3) break; // the line is too short; then stop                                                                                                   
    if (strncmp(p, "@RG\t", 4) == 0) {
      int c;
      khint_t k;
      if ((r = strstr(p, "\tID:")) != 0) {
	r += 4;
	for (s = r; *s != '\0' && *s != '\n' && *s != '\t'; ++s);
	c = *s; *s = '\0';
	k = kh_get(rg, h, r);
	*s = c;
	if (k != kh_end(h)) toprint = 1;
      }
    } else toprint = 1;
    if (toprint) {
      kputsn(p, q - p, &str); kputc('\n', &str);
    }
    p = q + 1;
  }
  *len = str.l;
  return str.s;
}

// ********* FOLLOWING BLOCK OF FUNCTIONS AND VARIABLES HAVE BEEN COPIED FROM bamc. and bgzf.c (and modified as needed) ******** //

static
void
report_error(BGZF* fp, const char* message) {
  fp->error = message;
}

static
int
deflate_block(BGZF* fp, int block_length)
{
  // Deflate the block in fp->uncompressed_block into fp->compressed_block.                                                                                         
  // Also adds an extra field that stores the compressed block length.                                                                                              

  bgzf_byte_t* buffer = fp->compressed_block;
  int buffer_size = fp->compressed_block_size;

  // Init gzip header                                                                                                                                               
  buffer[0] = GZIP_ID1;
  buffer[1] = GZIP_ID2;
  buffer[2] = CM_DEFLATE;
  buffer[3] = FLG_FEXTRA;
  buffer[4] = 0; // mtime                                                                                                                                           
  buffer[5] = 0;
  buffer[6] = 0;
  buffer[7] = 0;
  buffer[8] = 0;
  buffer[9] = OS_UNKNOWN;
  buffer[10] = BGZF_XLEN;
  buffer[11] = 0;
  buffer[12] = BGZF_ID1;
  buffer[13] = BGZF_ID2;
  buffer[14] = BGZF_LEN;
  buffer[15] = 0;
  buffer[16] = 0; // placeholder for block length
  buffer[17] = 0;

  // loop to retry for blocks that do not compress enough                                                                                                           
  int input_length = block_length;
  int compressed_length = 0;
  while (1) {
    z_stream zs;
    zs.zalloc = NULL;
    zs.zfree = NULL;
    zs.next_in = fp->uncompressed_block;
    zs.avail_in = input_length;
    zs.next_out = (void*)&buffer[BLOCK_HEADER_LENGTH];
    zs.avail_out = buffer_size - BLOCK_HEADER_LENGTH - BLOCK_FOOTER_LENGTH;

    int status = deflateInit2(&zs, fp->compress_level, Z_DEFLATED,
			      GZIP_WINDOW_BITS, Z_DEFAULT_MEM_LEVEL, Z_DEFAULT_STRATEGY);
    if (status != Z_OK) {
      report_error(fp, "deflate init failed");
      return -1;
    }
    status = deflate(&zs, Z_FINISH);
    if (status != Z_STREAM_END) {
      deflateEnd(&zs);
      if (status == Z_OK) {
	// Not enough space in buffer.                                                                                                                        
	// Can happen in the rare case the input doesn't compress enough.                                                                                     
	// Reduce the amount of input until it fits.                                                                                                          
	input_length -= 1024;
	if (input_length <= 0) {
	  // should never happen                                                                                                                            
	  report_error(fp, "input reduction failed");
	  return -1;
	}
	continue;
      }
      report_error(fp, "deflate failed");
      return -1;
    }
    status = deflateEnd(&zs);
    if (status != Z_OK) {
      report_error(fp, "deflate end failed");
      return -1;
    }
    compressed_length = zs.total_out;
    compressed_length += BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH;
    if (compressed_length > MAX_BLOCK_SIZE) {
      // should never happen                                                                                                                                    
      report_error(fp, "deflate overflow");
      return -1;
    }
    break;
  }

  packInt16((uint8_t*)&buffer[16], compressed_length-1);
  uint32_t crc = crc32(0L, NULL, 0L);
  crc = crc32(crc, fp->uncompressed_block, input_length);
  packInt32((uint8_t*)&buffer[compressed_length-8], crc);
  packInt32((uint8_t*)&buffer[compressed_length-4], input_length);

  int remaining = block_length - input_length;
  if (remaining > 0) {
    if (remaining > input_length) {
      // should never happen (check so we can use memcpy)                                                                                                       
      report_error(fp, "remainder too large");
      return -1;
    }
    memcpy(fp->uncompressed_block,
	   fp->uncompressed_block + input_length,
	   remaining);
  }
  fp->block_offset = remaining;
  return compressed_length;
}

static inline int bgzf_min(int x, int y)
{
  return (x < y) ? x : y;
}

static void swap_endian_data(const bam1_core_t *c, int data_len, uint8_t *data)
// copied one to one from bam.c
{
  uint8_t *s;
  uint32_t i, *cigar = (uint32_t*)(data + c->l_qname);
  s = data + c->n_cigar*4 + c->l_qname + c->l_qseq + (c->l_qseq + 1)/2;
  for (i = 0; i < c->n_cigar; ++i) bam_swap_endian_4p(&cigar[i]);
  while (s < data + data_len) {
    uint8_t type;
    s += 2; // skip key                                                                                                                                   
    type = toupper(*s); ++s; // skip type                                                                                                                 
    if (type == 'C' || type == 'A') ++s;
    else if (type == 'S') { bam_swap_endian_2p(s); s += 2; }
    else if (type == 'I' || type == 'F') { bam_swap_endian_4p(s); s += 4; }
    else if (type == 'D') { bam_swap_endian_8p(s); s += 8; }
    else if (type == 'Z' || type == 'H') { while (*s) ++s; ++s; }
    else if (type == 'B') {
      int32_t n, Bsize = bam_aux_type2size(*s);
      memcpy(&n, s + 1, 4);
      if (1 == Bsize) {
      } else if (2 == Bsize) {
	for (i = 0; i < n; i += 2)
	  for (i = 0; i < n; i += 2)
	    bam_swap_endian_2p(s + 5 + i);
      } else if (4 == Bsize) {
	for (i = 0; i < n; i += 4)
	  bam_swap_endian_4p(s + 5 + i);
      }
      bam_swap_endian_4p(s+1);
    }
  }
}

int bgzf_flush_try_x(BGZF *fp, int size,int *status)
{
  if (fp->block_offset + size > fp->uncompressed_block_size) {
    //return bgzf_flush(fp);
    *status = 1;
    return 0;
  }
  return -1;
}

inline int bam_write1_core_x(bamFile fp, const bam1_core_t *c, int data_len, uint8_t *data, bamFile fp2, int *status)
{
  uint32_t x[8], block_len = data_len + BAM_CORE_SIZE, y;
  int i;
  assert(BAM_CORE_SIZE == 32);
  x[0] = c->tid;
  x[1] = c->pos;
  x[2] = (uint32_t)c->bin<<16 | c->qual<<8 | c->l_qname;
  x[3] = (uint32_t)c->flag<<16 | c->n_cigar;
  x[4] = c->l_qseq;
  x[5] = c->mtid;
  x[6] = c->mpos;
  x[7] = c->isize;
  bgzf_flush_try_x(fp, 4 + block_len, status);
  if (bam_is_be) {
    for (i = 0; i < 8; ++i) bam_swap_endian_4p(x + i);
    y = block_len;
    // bam_write(fp, bam_swap_endian_4p(&y), 4);
    if(*status) bgzf_write_x(fp2, bam_swap_endian_4p(&y), 4, fp2, status); else bgzf_write_x(fp, bam_swap_endian_4p(&y), 4, fp2, status);
    swap_endian_data(c, data_len, data);
  } else // bam_write(fp, &block_len, 4); 
    if(*status) bgzf_write_x(fp2, &block_len, 4, fp2,status); else bgzf_write_x(fp, &block_len, 4, fp2,status);
  //bam_write(fp, x, BAM_CORE_SIZE);
  if (*status) bgzf_write_x(fp2, x, BAM_CORE_SIZE, fp2,status); else bgzf_write_x(fp, x, BAM_CORE_SIZE, fp2,status);
  //bam_write(fp, data, data_len);
  if(*status) bgzf_write_x(fp2, data, data_len, fp2,status); else bgzf_write_x(fp, data, data_len, fp2,status);
  if (bam_is_be) swap_endian_data(c, data_len, data);
  return 4 + block_len;
}

int bam_write1_x(bamFile fp, const bam1_t *b,bamFile fp2, int *status)
{
  return bam_write1_core_x(fp, &b->core, b->data_len, b->data, fp2, status);
}

int bgzf_write_x(BGZF* fp, const void* data, int length, BGZF* fp2, int *status)
{
  const bgzf_byte_t *input = data;
  int block_length, bytes_written;
  if (fp->open_mode != 'w') {
    report_error(fp, "file not open for writing");
    return -1;
  }
  
  if (fp->uncompressed_block == NULL)
    fp->uncompressed_block = malloc(fp->uncompressed_block_size);
  
  input = data;
  block_length = fp->uncompressed_block_size;
  bytes_written = 0;
  while (bytes_written < length) {
    int copy_length = bgzf_min(block_length - fp->block_offset, length - bytes_written);
    bgzf_byte_t* buffer = fp->uncompressed_block;
    memcpy(buffer + fp->block_offset, input, copy_length);
    fp->block_offset += copy_length;
    input += copy_length;
    bytes_written += copy_length;
    if (fp->block_offset == block_length) {
      //if (bgzf_flush(fp) != 0) {
      //	break;
      //}
      *status=1;
    }
  }
  return bytes_written;
}


// ********* END OF BLOCK OF FUNCTIONS AND VARIABLES COPIED/MODIFIED FROM bamc. and bgzf.c ******** //

typedef struct {
  int tid;
  BGZF *bf;
  char *ptr;
  size_t size;
} thread_aux_t;

static void *flush_worker_thread(void *data) {

  thread_aux_t *d = (thread_aux_t*)data;
  bgzf_flush(d->bf);

  return 0;
}

int flush_worker(BGZF *bf) {
  bgzf_flush(bf);
  fflush(bf->file);
  fclose(bf->file);

  return 0;
}

int flush(BGZF *bf[], int bflen, int start, int end, int elem, FILE *fh_out) {
  // bf - array of BGZF structs
  // bflen - number of elements in bf array
  // start - first position in array to be processed/flushed
  // end - last position in array to be processed/flushed
  // elem - nr of elements from start to end (exclusive) to be processed.
  // The array will be processed from 'start' position one by one towards the bflen of the array.
  // Tn case bflen is reached, processing will continue from 0 towards bflen.
  // The last position to be processed is 'end-1'.

  int i;
  char *ptr[bflen];
  size_t size[bflen];

  // assign pthread related variables
  pthread_attr_t attr;
  thread_aux_t *data;
  pthread_t tid[elem];
  int tidnr = 0;

  // create array of thread_aux_t used to hold the data of threads
  data = (thread_aux_t*)calloc(elem, sizeof(thread_aux_t));

  // initialise pthreads
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  // flush in threads
  i = start;
  while(1) {
    if(i == bflen) i = 0;
    if(i == end ) break;

    // perform bgzf_flush (includes data compression we are interested in) in threads
    data[tidnr].tid = tidnr;
    data[tidnr].bf = bf[i];
    data[tidnr].bf->file = open_memstream(&data[tidnr].ptr, &data[tidnr].size);
    pthread_create(&tid[tidnr],&attr,flush_worker_thread,data+tidnr);
    tidnr++;
    
    // An example of doing the same as above without threads
    // flush_worker(bf[i]);
    // fwrite(ptr[i], 1, size[i], fh_out);
    // free(ptr[i]);
    // size[i] = 0;
    
    i++;
  }

  // join the threads
  i = start;
  tidnr = 0;
  while(1) {
    if(i == bflen) i = 0;
    if(i == end ) break;
    
    pthread_join(tid[tidnr], 0);
    fclose(data[tidnr].bf->file);
    fwrite(data[tidnr].ptr, 1, data[tidnr].size, fh_out);
    free(data[tidnr].ptr);

    tidnr++;
    i++;
  }
  free(data);

  return 0;
}

int sam2bam(char *fn_in, char *fn_out, int nthreads) {

  // fn_in  - name of the input sam file
  // fn_out - name of the output bam file. If '-' the file will be printed to stdout.

  int ret = 0, compress_level = -1, is_count = 0, is_header=1;
  int count = 0;
  samfile_t *in = 0, *out = 0;
  char in_mode[5], out_mode[5];
  FILE *fh_out;

  strcpy(in_mode, "r");
  // set output mode to bam
  strcpy(out_mode, "w");
  strcat(out_mode, "b");
  if (is_header) strcat(out_mode, "h");

  // set bam compression level to default compression
  compress_level = 0;

  // open input file and read header.
  if ((in = samopen(fn_in, in_mode, 0)) == 0) {
    fprintf(stderr, "[main_samview] fail to open \"%s\" for reading.\n", fn_in);
    ret = 1;
  }
  
  // check if header OK
  if (in->header == 0) {
    fprintf(stderr, "[main_samview] fail to read the header from \"%s\".\n", fn_in);
    ret = 1;
  }

  // copied from samtools. No idea what this block does. HL hints in his comments that he does not either.
  if (g_rghash) { // FIXME: I do not know what "bam_header_t::n_text" is for...
    char *tmp;
    int l;
    tmp = drop_rg(in->header->text, g_rghash, &l);
    free(in->header->text);
    in->header->text = tmp;
    in->header->l_text = l;
  }

  // open output
  if(strcmp(fn_out,"-")) fh_out = fopen(fn_out,"w");
  else fh_out = stdout;  
  if(fh_out == 0) {
    fprintf(stderr, "[main_samview] fail to open \"%s\" for writing.\n", fn_out? fn_out : "standard output");
    ret = 1;
  }

  // return in case we've failed
  if(ret) {
    return ret;
  }

  int inthreads = nthreads+1; // internal number of threads is used to assign and operate the number of BGZF structs holding pieces of data during compression.
  int cthread = 0;
  int cthread2 = 1;
  char *ptr;
  size_t size;
  BGZF* bf[inthreads];
  int i;
  int status = 0;
  int follower = 0;
  
  // install BGZF struct. TODO: The way of how the bgzf structs are installed via having to open a file has to be changed.
  for(i=0;i<inthreads;i++) {
    fprintf(stderr,"Open\n");
    bf[i] = bgzf_open("tmp.hohoho", "w");
    //fprintf(stderr,"Openened\n");
    fclose(bf[i]->file);
    // fprintf(stderr,"Redirected\n");
  }
  
  fprintf(stderr,"[sam2bam] writing sam header\n");
  // write bam header
  // TODO: rather than installing another bgzf, could use the first one from the bgzf array instead.
  BGZF* bfh;
  char *ptrh;
  size_t sizeh;
  bfh = bgzf_open("tmp.hohoho", "w");
  fclose(bfh->file);
  bfh->file = open_memstream(&ptrh, &sizeh);
  bam_header_write(bfh,in->header);  
  fclose(bfh->file);
  fwrite(ptrh, 1, sizeh, fh_out);
  free(ptrh);

  int nrypes = 0;

  // in case not failed so far, process the rest of the file one read at the time.
  if(!ret) {
    bam1_t *b = bam_init1();
    int r;
    while ((r = samread(in, b)) >= 0) { // read one alignment from `in'                                                 
      if (!__g_skip_aln(in->header, b)) {
	if (!is_count) {
	  
	  // write sequence but without doing the compression quite yet
	  bam_write1_x(bf[cthread],b,bf[cthread2],&status);
	  // check if block is ready for compression and flushing on disk
	  if(status == 1) {	    
	    nrypes++;

	    // shift the active BGZF to next one
	    cthread = cthread2;
	    if((cthread2+1) == inthreads) cthread2=0; else cthread2++;
	    status = 0;
	    
	    // when enough blocks ready for compression/flushing, process these in parallel
	    if(nrypes == (inthreads-1)) {
	      flush(bf,inthreads,cthread+1,cthread, nrypes, fh_out);
	      nrypes = 0;
	    }

	  } // end if status == 1

	}
	count++;
      }
    }
    if (r < -1) {
      fprintf(stderr, "[main_samview] truncated file.\n");
      ret = 1;
    }
    bam_destroy1(b);
  }

  // compress/flush the remaining blocks
  if(nrypes != 0) {
    i = cthread - nrypes;
    if(i < 0) i = inthreads + i;  
    flush(bf, inthreads, i, cthread, nrypes, fh_out);
  }
  // compress/flush (done automatically on bgzf_close) the last remaining block
  bf[cthread]->file = open_memstream(&ptr, &size);
  bgzf_close(bf[cthread]);
  fwrite(ptr, 1, size, fh_out);
  free(ptr);
  size = 0;  

  // close output file
  fclose(fh_out);

  fprintf(stderr,"[sam2bam] %d reads processed.\n",count);
  
  // do the aftermath
  if (is_count && ret == 0) {
    printf("%d\n", count);
  }

  // close files, free and return
  free(g_rg);
  if (g_bed) bed_destroy(g_bed);
  if (g_rghash) {
    khint_t k;
    for (k = 0; k < kh_end(g_rghash); ++k)
      if (kh_exist(g_rghash, k)) free((char*)kh_key(g_rghash, k));
    kh_destroy(rg, g_rghash);
  }
  samclose(in);
  if (!is_count)
    samclose(out);
  return ret;
}

int main(int argc, char *argv[]){
  int nthreads = 1;
  int c;

  while ((c = getopt(argc, argv, "t:")) >= 0) {
    switch (c) {
    case 't':
      nthreads = atoi(optarg);
      break;
    default: return 1;
    }
  }

  if ((optind + 2) != argc) {
    fprintf(stderr,"\nUsage:   sam2bam [-t] <filename.sam> <filename.bam>\n");
    fprintf(stderr,"\nOptions: -t number of threads (CPU cores) to use (default = 1)\n\n");
    return 1;
  }

  return sam2bam(argv[optind],argv[optind+1], nthreads);
}

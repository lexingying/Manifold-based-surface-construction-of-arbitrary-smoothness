/****************************************************************/
/*    NAME: David Powell                                        */
/*    ACCT: dep                                                 */
/*    FILE: IImage.C                                            */
/****************************************************************/

#include "IImage.hpp"

#include <iostream>
//using namespace std;

#include <png.h>

// for: #define offsetof(TYPE, MEMBER) ((size_t) &((TYPE *)0)->MEMBER)
// we use this to determine the pixel format on the fly
#include <stddef.h>
#include <assert.h>

/* modified by LEXING YING for the fluid-shell code (voxelization) */

using namespace std;

// ======================================================================
//
//  Function: IImage::IImage
//     Takes: 
//   Returns: nothing
//   Effects: 
//
// ======================================================================

IImage::IImage() : m_width(-1), m_height(-1), m_data(NULL)
{
  // assert( sizeof(IPixel) == 4 ) ;
  if( sizeof(IPixel) != 4 ) {
    fprintf( stderr, "sizeof(IPixel) != 4.  Struct padding?\n" ) ;
    exit(1) ;
  }
}

IImage::IImage( const IImage& copy )
{
  m_width = copy.m_width ;
  m_height = copy.m_height ;
  m_data = (IPixel*) malloc( m_width * m_height * sizeof(IPixel) ) ;
  memcpy( m_data, copy.m_data, m_width * m_height * sizeof(IPixel) ) ;
}


// ======================================================================
//
//  Function: IImage::~IImage
//     Takes: void
//   Returns: nothing
//   Effects: Destroys the window, cleaning up as necessary
//
// ======================================================================

IImage::~IImage()
{
  if(m_data) free(m_data);
  m_data = NULL;
}


// We use out own reading/writing functions because libpng may have been compiled using a
// different compiler & libc.  This could make the FILE*'s incompatible.
static void user_write_data(png_structp png_ptr, png_bytep data, png_size_t length);
static void user_write_data(png_structp png_ptr, png_bytep data, png_size_t length)
{
  voidp write_io_ptr = png_get_io_ptr(png_ptr) ;
  fwrite( (unsigned char*) data, length, 1, (FILE*) write_io_ptr ) ;
}

static void user_flush_data(png_structp png_ptr);
static void user_flush_data(png_structp png_ptr)
{
  voidp write_io_ptr = png_get_io_ptr(png_ptr) ;
  fflush( (FILE*) write_io_ptr ) ;
}

static void user_read_data(png_structp png_ptr, png_bytep data, png_size_t length) ;
static void user_read_data(png_structp png_ptr, png_bytep data, png_size_t length)
{
  voidp read_io_ptr = png_get_io_ptr(png_ptr) ;
  fread( (unsigned char*) data, length, 1, (FILE*) read_io_ptr ) ;
}


bool IImage::save(const char* fname)
{
  FILE* fp = NULL ;
  bool rval = true ;
  png_structp png_ptr = NULL ;
  png_infop info_ptr = NULL ;
  
  // Some error strings
  char open[] = "Error opening %s\n" ;
  char problem[] = "libpng encountered a problem writing %s)\n" ;
  
  // Open the file for reading in binary mode.
  fp = fopen( fname, "wb" ) ;
  if( !fp ) {
    fprintf( stderr, open, fname ) ;
    rval = false ;
    goto IImage_save_cleanup ;
  }
  
  // Allocate the png structs.
  png_ptr = png_create_write_struct( PNG_LIBPNG_VER_STRING, NULL, NULL, NULL ) ;
  if( !png_ptr ) {
    fprintf( stderr, problem, fname ) ;
    rval = false ;
    goto IImage_save_cleanup ;
  }
  info_ptr = png_create_info_struct( png_ptr ) ;
  if( !info_ptr ) {
    fprintf( stderr, problem, fname ) ;
    rval = false ;
    goto IImage_save_cleanup ;
    }
  
  // Set up the png error routine.
  if( setjmp(png_jmpbuf(png_ptr)) ) {
    fprintf( stderr, problem, fname ) ;
    rval = false ;
    goto IImage_save_cleanup ;
    }
  
  // Give libpng the FILE*.
  // png_init_io( png_ptr, fp ) ;
  // or
  // use our own write callback
  png_set_write_fn( png_ptr, fp, (png_rw_ptr) user_write_data, user_flush_data ) ;
  
  // We'll use the low-level interface since the high-level interface won't handle
  // png_set_filler() which we need to tell libpng to strip out the A from our
  // 4-byte pixels.
  
  // First we set and write the info struct.
  png_set_IHDR( png_ptr, info_ptr,
    m_width, m_height,
    8, PNG_COLOR_TYPE_RGB,
    PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT ) ;
  png_write_info( png_ptr, info_ptr ) ;
  
  // Then we set up transforms.
  // 1. tell libpng to strip out the filler byte in our 4-byte pixels
  // if IPixel::a comes after any other member (b,g,r), we stip AFTER
  if( offsetof( IPixel, a ) > offsetof( IPixel, b ) ) {
    png_set_filler( png_ptr, 0, PNG_FILLER_AFTER ) ;
    // printf("alpha after\n");
  } else {
    png_set_filler( png_ptr, 0, PNG_FILLER_BEFORE ) ;
    // printf("alpha before\n");
  }
  // 2. tell libpng how our color triples are stored (b < r or vice versa)
  if( offsetof( IPixel, b ) < offsetof( IPixel, r ) ) {
    png_set_bgr( png_ptr ) ;
    // printf("bgr\n") ;
  }
  // else { printf("rgb\n"); }
  
  // Finally we create a row_pointers[] pointing into our data* and write the png out to the FILE*.
  {
    // 1. allocate row pointers array
    png_bytep* row_pointers = (png_bytep*) png_malloc( png_ptr, m_height * sizeof(png_bytep) ) ;
    // 2. point row pointers into m_data
    for( int i = 0 ; i < m_height ; ++i ) {
      row_pointers[i] = (png_bytep) (m_data + i*m_width) ;
    }
    // 3. write the image data
    png_write_image( png_ptr, row_pointers ) ;
    // 4. free row pointers array
    png_free( png_ptr, row_pointers ) ;
  }
  
  // Write out end info.  We're done.  Fall through to cleanup.
  png_write_end( png_ptr, NULL ) ;
  
IImage_save_cleanup:
  png_destroy_write_struct( png_ptr ? &png_ptr : NULL, info_ptr ? &info_ptr : NULL ) ;
  if( fp ) fclose( fp ) ;
  
  return rval ;
}

bool IImage::load(const char* fname)
{
  // TODO use user_read_data instead of libpng's... safer (took it out
  // for sake of simplicity)
  FILE *pfile = fopen(fname, "rb");
  if (NULL==pfile) {
    return false;
  }
  png_structp png_ptr = png_create_read_struct
    (PNG_LIBPNG_VER_STRING, (png_voidp)NULL,
     NULL, NULL);
  if (NULL==png_ptr) {
    return false;
  }
  png_infop info_ptr = png_create_info_struct(png_ptr);
  if (NULL==info_ptr) {
    png_destroy_read_struct(&png_ptr,
        (png_infopp)NULL, (png_infopp)NULL);
    return false;
  }

  png_infop end_info = png_create_info_struct(png_ptr);
  if (NULL==end_info) {
    png_destroy_read_struct(&png_ptr,
        &info_ptr, (png_infopp)NULL);
    return false;
  }

  if (setjmp(png_jmpbuf(png_ptr)))
  {
    png_destroy_read_struct(&png_ptr, &info_ptr,
        &end_info);
    fclose(pfile);
    return false;
  }


  /* png_ptr, info_ptr, and end_info are initialized at this point */

  // In case libpng was compiled with a different stdio:
  // png_init_io(png_ptr, pfile);
  png_set_read_fn( png_ptr, pfile, (png_rw_ptr) user_read_data ) ;

  /* this is where we would call png_set_sig_bytes() */

  /* this is where we would call read_row_callback() */

  {
    int png_transforms = 
      PNG_TRANSFORM_STRIP_16 |
      PNG_TRANSFORM_STRIP_ALPHA |
      PNG_TRANSFORM_PACKING |
      PNG_TRANSFORM_PACKSWAP |
      PNG_TRANSFORM_EXPAND |
      PNG_TRANSFORM_SHIFT;
    png_read_png(png_ptr, info_ptr, png_transforms, NULL);

    png_uint_32 width = png_get_image_width(png_ptr, info_ptr);
    png_uint_32 hei = png_get_image_height(png_ptr, info_ptr);
    png_uint_32 bytes = png_get_rowbytes(png_ptr, info_ptr);
    //png_uint_32 bitdepth = png_get_bit_depth(png_ptr, info_ptr);
    //png_uint_32 coltype = png_get_color_type(png_ptr, info_ptr);
    png_byte channels = png_get_channels(png_ptr, info_ptr);

    png_bytepp row_pointers;
    row_pointers = png_get_rows(png_ptr, info_ptr);
    
    assert(NULL!=row_pointers);
#if 0
    printf("width %d, hei %d, rowbytes %d, bitd %d, "
      "colt %d, channels %d\n",
      (int)width,(int)hei,(int)bytes,(int)bitdepth,
      (int)coltype,(int)channels);
#endif

    resize(width,hei);

    // TODO this is so gross.  I should be in jail.
    unsigned int rowc;
    unsigned int hackcount;
    const unsigned char *rp;
    for (rowc=0; rowc<hei; rowc++) {
      if (1 == channels) {
        // if mono, must expand to fill all four bytes
        // for the time being
        assert(bytes/width == 1);
        IPixel *mediary = m_data + width*rowc;
        rp = row_pointers[rowc];
        // expand each sample into BG&R, then add full
        // alpha:
        // (again, I would like to see a transform for
        // this, but there isn't.  The transforms do
        // not change the number of channels, it seems)
        for (hackcount=0; hackcount<bytes; hackcount++) {
          mediary[hackcount].r = rp[hackcount];
          mediary[hackcount].g = rp[hackcount];
          mediary[hackcount].b = rp[hackcount];
          mediary[hackcount].a = 255;
        }
      } else if (3 == channels) {
        // if color, must expand to fill all four bytes
        // for the time being
        assert(bytes/width == 3);
        IPixel *mediary = m_data + width*rowc;
        rp = row_pointers[rowc];
        int idx = 0;
        // the png image has no alpha.  we do.  add a
        // byte.
        // (there should be a transform filter
        // for this above, but there isn't)
        // because of the transform, we know that
        // libpng is giving us consecutive "rgb" triples.
        for (hackcount=0; hackcount<width; hackcount++) {
          mediary[hackcount].r = rp[idx++];
          mediary[hackcount].g = rp[idx++];
          mediary[hackcount].b = rp[idx++];
          mediary[hackcount].a = 255;
        }
      } else {
        fprintf(stderr, "unsupported channel depth in"
            " png \"%s\".  Exiting.\n",
            fname);
      }
    }
  }

  png_destroy_read_struct(&png_ptr, &info_ptr, &end_info);
  fclose(pfile);
  return true;
}




// ======================================================================
//
//  Function: IImage::data
//     Takes: void
//   Returns: usigned char* - a pointer to the beginning of the raw image
//            data
//   Effects: obtains the aforementioned pointer
//
// ======================================================================

IPixel* IImage::data()
{
  return m_data;
}





// ======================================================================
//
//  Function: IImage::width
//     Takes: void
//   Returns: int - the width of the image
//   Effects: obtains the width
//
// ======================================================================

int IImage::width() const
{
  return m_width;
}





// ======================================================================
//
//  Function: IImage::height
//     Takes: void
//   Returns: int - the height of the image
//   Effects: obtains the height
//
// ======================================================================

int IImage::height() const
{
  return m_height;
}





// ======================================================================
//
//  Function: IImage::resize
//     Takes: width,height - the new width and height
//   Returns: void
//   Effects: Destroys the existing image and creates a new one with the
//            specified size
//
// ======================================================================

void IImage::resize( int width, int height )
{
  if( m_data ) free( m_data ) ;
  m_width = width ;
  m_height = height ;
  m_data = (IPixel*) malloc( m_width * m_height * sizeof(IPixel) ) ;
  memset( m_data, 0, m_width * m_height * sizeof(IPixel) ) ;
}





// ======================================================================
//
//  Function: IImage::putPixel
//     Takes: x,y - the pixel to modify
//            r,g,b - the color to set the pixel from
//   Returns: void
//   Effects: Changes the color of the specified pixel
//
// ======================================================================

void IImage::putPixel(int x, int y, int r, int g, int b)
{
  if(x<0||y<0||x>=m_width||y>=m_height) return;
  IPixel* pixel = m_data + (y*m_width+x) ;
  pixel->r = r;
  pixel->g = g;
  pixel->b = b;
}





// ======================================================================
//
//  Function: IImage::putRed
//     Takes: x,y - the pixel to modify
//            value - the new red value
//   Returns: void
//   Effects: changes the red value of the specified pixel
//
// ======================================================================

void IImage::putRed(int x, int y, int value)
{
  if(x<0||y<0||x>=m_width||y>=m_height) return;
  m_data[ y*m_width+x ].r = value;
}


IPixel*
IImage::getPixel(int x, int y)
{
  assert(!(x<0||y<0||x>=m_width||y>=m_height));
  return (m_data+(y*m_width+x));
}

const IPixel*
IImage::getPixel(int x, int y) const {
  assert(!(x<0||y<0||x>=m_width||y>=m_height));
  return (m_data+(y*m_width+x));
}



// ======================================================================
//
//  Function: IImage::putGreen
//     Takes: x,y - the pixel to modify
//            value - the new green value
//   Returns: void
//   Effects: changes the green value of the specified pixel
//
// ======================================================================

void IImage::putGreen(int x, int y, int value)
{
  if(x<0||y<0||x>=m_width||y>=m_height) return;
  m_data[ y*m_width+x ].g = value;
}





// ======================================================================
//
//  Function: IImage::putBlue
//     Takes: x,y - the pixel to modify
//            value - the new blue value
//   Returns: void
//   Effects: changes the blue value of the specified pixel
//
// ======================================================================

void IImage::putBlue(int x, int y, int value)
{
  if(x<0||y<0||x>=m_width||y>=m_height) return;
  m_data[ y*m_width+x ].b = value;
}





// ======================================================================
//
//  Function: IImage::getRed
//     Takes: x,y - the pixel to access
//   Returns: int - the red value of the pixel
//   Effects: 
//
// ======================================================================

int IImage::getRed(int x, int y)
{
  if(x<0||y<0||x>=m_width||y>=m_height) return -1;
  return m_data[ y*m_width+x ].r ;
}





// ======================================================================
//
//  Function: IImage::getGreen
//     Takes: x,y - the pixel to access
//   Returns: int - the green value of the pixel
//   Effects: 
//
// ======================================================================

int IImage::getGreen(int x, int y)
{
  if(x<0||y<0||x>=m_width||y>=m_height) return -1;
  return m_data[ y*m_width+x ].g ;
}





// ======================================================================
//
//  Function: IImage::getBlue
//     Takes: x,y - the pixel to access
//   Returns: int - the blue value of the pixel
//   Effects: 
//
// ======================================================================

int IImage::getBlue(int x, int y)
{
  if(x<0||y<0||x>=m_width||y>=m_height) return -1;
  return m_data[ y*m_width+x ].b ;
}






#if 0
// YOTAM'S OLD LOAD(), trying mine in order to keep it simple:
  FILE* fp = NULL ;
  bool rval = true;
  png_structp png_ptr = NULL ;
  png_infop info_ptr = NULL ;
  
  // Some error strings
  char open[] = "Error opening %s\n" ;
  char invalid[] = "Invalid png file: %s\n" ;
  char liberr[] = "libpng encountered a problem (not related to file)\n" ;
  char maybe[] = "libpng encountered a problem (may be related to file: %s)\n" ;
  
  // for checking the png header
  const size_t PNG_BYTES_TO_CHECK = 4 ; // example.c uses 4
  png_byte header[ PNG_BYTES_TO_CHECK ] ;
  
  // Open the file for reading in binary mode.
  fp = fopen( fname, "rb" ) ;
  if( !fp ) {
    fprintf( stderr, open, fname ) ;
    rval = false;
    goto IImage_load_cleanup ;
  }
  
  // Check some bytes at the beginning of the file to make sure it's a png.
  if( PNG_BYTES_TO_CHECK != fread( header, 1, PNG_BYTES_TO_CHECK, fp ) ) {
    fprintf( stderr, invalid, fname ) ;
    rval = false;
    goto IImage_load_cleanup ;
  }
  if( png_sig_cmp( header, 0, PNG_BYTES_TO_CHECK ) ) {
    fprintf( stderr, invalid, fname ) ;
    rval = false;
    goto IImage_load_cleanup ;
  }
  
  // Since it looks like we have a good png file, allocate the png structs.
  png_ptr = png_create_read_struct( PNG_LIBPNG_VER_STRING, NULL, NULL, NULL ) ;
  if( !png_ptr ) {
    fprintf( stderr, liberr ) ;
    rval = false;
    goto IImage_load_cleanup ;
  }
  info_ptr = png_create_info_struct( png_ptr ) ;
  if( !info_ptr ) {
    fprintf( stderr, liberr ) ;
    rval = false;
    goto IImage_load_cleanup ;
    }
  
  // Set up the png error routine.
  if( setjmp(png_jmpbuf(png_ptr)) ) {
    fprintf( stderr, maybe, fname ) ;
    rval = false;
    goto IImage_load_cleanup ;
    }
  
  // Give libpng the FILE*, tell it how many bytes we read.
  // png_init_io( png_ptr, fp ) ;
  // or
  // use our own read callback
  png_set_read_fn( png_ptr, fp, (png_rw_ptr) user_read_data ) ;
  
  png_set_sig_bytes( png_ptr, PNG_BYTES_TO_CHECK ) ;
  
  // We'll use the low-level interface since the high-level interface won't handle
  // png_set_filler() which we need to guarantee there'll be a filler "A" in our
  // 4-byte ARGB pixels.
  // Really the low-level interface isn't more complicated than the high-level interface.
  // To choose transform flags we have to query the png_info struct.
  // Instead of OR'ing in another transform flag (high-level interface), we call a set
  // method (low-level interface).
  
  // First we read the info struct.
  png_read_info( png_ptr, info_ptr ) ;
  
  // Now we set up transforms.
  // 1. convert gray and paletted to rgb (this guarantees 8 or 16 bit depths (rgb must be 8 or 16))
  // 2. convert 16 bit depths to 8
  // 3. super-fanciness: we set gamma and add a background color for images
  //    with an alpha channel (currently off)
  // 4. if we don't have alpha, add a filler before, otherwise prepend alpha to the color triple
  // 5. If we were on big-endian, ask for BGR (instead of RGB) here
  // 6. ask libpng to deinterlace
  {
    
    // 1
    png_byte color_type = png_get_color_type( png_ptr, info_ptr ) ;
    if( color_type == PNG_COLOR_TYPE_PALETTE )
      png_set_expand( png_ptr ) ;
    else if( color_type == PNG_COLOR_TYPE_GRAY || color_type == PNG_COLOR_TYPE_GRAY_ALPHA )
      png_set_gray_to_rgb( png_ptr ) ;
    
    // 2
    // png_read_update_info( png_ptr, info_ptr ) ;
    png_byte depth = png_get_bit_depth( png_ptr, info_ptr ) ;
    // assert( depth >= 8 ) ; // since we call png_read_update_info()
    
    if( 8 < depth )
      png_set_strip_16( png_ptr ) ;
    
#define FANCY_SHMANCY 0
#if FANCY_SHMANCY
    // 3 (code taken from png's example.c)
    // background
    // NOTE: Setting background must occur before the alpha/filler code.
    // NOTE: If background-setting happens afterwards, it will remove
    // NOTE: the alpha channel and less 1 byte from each pixel.
    {
      // png_color_16 default_background = { 0, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF } ; // white
      png_color_16 default_background = { 0, 0,0,0,0 } ; // black
      png_color_16p image_background ;
      
      if( png_get_bKGD( png_ptr, info_ptr, &image_background ) )
        png_set_background( png_ptr, image_background, PNG_BACKGROUND_GAMMA_FILE, 1, 1.0 ) ;
      else
            png_set_background( png_ptr, &default_background, PNG_BACKGROUND_GAMMA_SCREEN, 0, 1.0 ) ;
    }
    
    // gamma
    {
      double screen_gamma ;
      {
        char* gamma_str ;
        if( (gamma_str = getenv("SCREEN_GAMMA")) != NULL )
          screen_gamma = atof( gamma_str ) ;
        else {
          // TODO: this is a bad assumption
          #if defined(__linux__) || defined(WIN32)
          screen_gamma = 2.2 ;
          #else
          screen_gamma = 1.7 ;
          #endif
        }
      }
      
      int intent ;
      if( png_get_sRGB( png_ptr, info_ptr, &intent ) )
        png_set_gamma( png_ptr, screen_gamma, 0.45455 ) ;
      else {
        double image_gamma ;
        if( png_get_gAMA( png_ptr, info_ptr, &image_gamma ) )
          png_set_gamma( png_ptr, screen_gamma, image_gamma ) ;
        else
          png_set_gamma( png_ptr, screen_gamma, 0.45455 ) ;
      }
    }
#endif
#undef FANCY_SHMANCY
    
    // NOTE: These next two affect the layout of the channels in the pixel
    // NOTE: Steps 4 and 5 turn the pixel into ABGR, but this can easily
    // NOTE: be changed to RGBA by
    // NOTE: { PNG_FILLER_BEORE -> PNG_FILLER_AFTER / not calling set_swap_alpha }
    // NOTE: and not calling set_bgr
    // FOLLOWUP: in fact, on linux we want ARGB, so we don't call set_bgr
    
    // 4
    if( color_type != PNG_COLOR_TYPE_GRAY_ALPHA && color_type != PNG_COLOR_TYPE_RGB_ALPHA ) {
    /*
    // since we call png_read_update_info(), we can use channels instead of the original pixel type
    png_read_update_info( png_ptr, info_ptr ) ;
    assert( png_get_color_type( png_ptr, info_ptr ) == PNG_COLOR_TYPE_RGB ||
      png_get_color_type( png_ptr, info_ptr ) == PNG_COLOR_TYPE_RGB_ALPHA ) ;
    png_byte channels = png_get_channels( png_ptr, info_ptr ) ;
    assert( channels == 3 || channels == 4 ) ;
    if( channels < 4 ) {
    */
      png_set_filler( png_ptr, 0xFF, PNG_FILLER_BEFORE ) ;
    } else {
      png_set_swap_alpha( png_ptr ) ;
    }
    
    // 5
    // png_set_bgr( png_ptr ) ;
    // 6
    png_set_interlace_handling( png_ptr ) ;
    
  }
  
  // We're almost ready to copy over the png data.
  // First we must resize our data* and set our width & height.
  {
    png_uint_32 width = png_get_image_width( png_ptr, info_ptr ) ;
    png_uint_32 height = png_get_image_height( png_ptr, info_ptr ) ;
    
    // fprintf( stderr, "width: %d, height: %d\n", (int) width, (int) height ) ;
    
    resize( width, height ) ;
  }
  
  // Now we can create a rows[] pointing into our data* and read the png into our buffer.
  {
    // fprintf( stderr, "width: %d, height: %d\n", (int) m_width, (int) m_height ) ;
    
    // 1. allocate row pointers array
    png_bytep* row_pointers = (png_bytep*) png_malloc( png_ptr, m_height * sizeof(png_bytep) ) ;
    // 2. point row pointers into m_data
    for( int i = 0 ; i < m_height ; ++i ) {
      row_pointers[i] = (png_bytep) (m_data + i*m_width ) ;
    }
    // 3. read the image data
    png_read_image( png_ptr, row_pointers ) ;
    // 4. free row pointers array
    png_free( png_ptr, row_pointers ) ;
  }
  
  // Read the end info.  We're done.  Fall through to cleanup.
  png_read_end( png_ptr, NULL ) ;
  
IImage_load_cleanup:
  // due to (what looks like) a bug in libpng-1.2.4, we can't pass NULL for the png_ptr arg
  if( png_ptr ) png_destroy_read_struct( &png_ptr, info_ptr ? &info_ptr : NULL, NULL ) ;
  
  if( fp ) fclose( fp ) ;

#if 0
  for (int c1=0; c1<m_height*m_width*4; c1++) {
    printf("%d ",(int)m_data[c1]);
  }
#endif
  
  return rval ;

#endif

/* modified by LEXING YING for the fluid-shell code (voxelization) */


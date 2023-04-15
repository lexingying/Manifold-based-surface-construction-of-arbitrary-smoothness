/****************************************************************/
/*    NAME: David Powell                                        */
/*    ACCT: dep                                                 */
/*    FILE: IImage.H                                            */
/****************************************************************/


/* modified by LEXING YING */
#ifndef _IIMAGE_HPP_
#define _IIMAGE_HPP_

#include "common/commoninc.hpp"

//#ifndef _IImage_H_
//#define _IImage_H_

// file loading/saving automatically picks up changes to this struct.
struct IPixel {
  unsigned char b ;
  unsigned char g ;
  unsigned char r ;
  unsigned char a ;
};

// ======================================================================
//
//  Class: IImage
//
// ======================================================================

class IImage {

protected:
  int m_width, m_height;
  IPixel* m_data; // The raw image data

public:
  IImage();
  IImage( const IImage& );
  virtual ~IImage();

  bool save(const char* fname);
  bool load(const char* fname);

  IPixel* data();

  int width() const;
  int height() const;
  virtual void resize(int width, int height);

  void putPixel(int x, int y, int r, int g, int b); 

  void putRed(int x, int y, int value);
  void putGreen(int x, int y, int value);
  void putBlue(int x, int y, int value);

  int getRed(int x, int y);
  int getGreen(int x, int y);
  int getBlue(int x, int y);

  IPixel* getPixel(int x, int y);
  const IPixel* getPixel(int x, int y) const;

};
/* modified by LEXING YING for the fluid-shell code (voxelization) */
//}

#endif


/*
 * 8-bit Gray Scale BMP loader and saver.
 *
 * This is an adaptation of bitmap.h to fit 8-bit gray .bmp files.
 * 
 * Adapted by: Ivo Knottnerus
 */

#ifndef BITMAP_H
#define BITMAP_H

#include <iostream>
#include <fstream>
#include <string>

#ifndef __LITTLE_ENDIAN__
	#ifndef __BIG_ENDIAN__
		#define __LITTLE_ENDIAN__
	#endif
#endif

#ifdef __LITTLE_ENDIAN__
	#define BITMAP_SIGNATURE 0x4d42
#else
	#define BITMAP_SIGNATURE 0x424d
#endif

#if defined(_MSC_VER) || defined(__INTEL_COMPILER)
	typedef unsigned __int32 uint32_t;
	typedef unsigned __int16 uint16_t;
	typedef unsigned __int8 uint8_t;
	typedef __int32 int32_t;
#elif defined(__GNUC__) || defined(__CYGWIN__) || defined(__MWERKS__) || defined(__WATCOMC__) || defined(__PGI) || defined(__LCC__)
	#include <stdint.h>
#else
	typedef unsigned int uint32_t;
	typedef unsigned short int uint16_t;
	typedef unsigned char uint8_t;
	typedef int int32_t;
#endif

#pragma pack(push, 1)

typedef struct _BITMAP_FILEHEADER {
	uint16_t Signature;
	uint32_t Size;
	uint32_t Reserved;
	uint32_t BitsOffset;
} BITMAP_FILEHEADER;

#define BITMAP_FILEHEADER_SIZE 14

typedef struct _BITMAP_HEADER {
	uint32_t HeaderSize;
	int32_t Width;
	int32_t Height;
	uint16_t Planes;
	uint16_t BitCount;
	uint32_t Compression;
	uint32_t SizeImage;
	int32_t PelsPerMeterX;
	int32_t PelsPerMeterY;
	uint32_t ClrUsed;
	uint32_t ClrImportant;
	uint32_t RedMask;
	uint32_t GreenMask;
	uint32_t BlueMask;
	uint32_t AlphaMask;
	uint32_t CsType;
	uint32_t Endpoints[9]; // see http://msdn2.microsoft.com/en-us/library/ms536569.aspx
	uint32_t GammaRed;
	uint32_t GammaGreen;
	uint32_t GammaBlue;
} BITMAP_HEADER;

typedef struct _BITMAP_HEADER_SHORT {
	uint32_t HeaderSize;
	int32_t Width;
	int32_t Height;
	uint16_t Planes;
	uint16_t BitCount;
	uint32_t Compression;
	uint32_t SizeImage;
	int32_t PelsPerMeterX;
	int32_t PelsPerMeterY;
	uint32_t ClrUsed;
	uint32_t ClrImportant;
} BITMAP_HEADER_SHORT;

#pragma pack(pop)

class CBitmap {
private:
	BITMAP_FILEHEADER m_BitmapFileHeader;
	BITMAP_HEADER_SHORT m_BitmapHeader;
	uint8_t *m_BitmapData;
	unsigned int m_BitmapSize;

public:
	
	CBitmap() : m_BitmapData(0), m_BitmapSize(0)  {
		Dispose();
	}
	
	CBitmap(const char* Filename) : m_BitmapData(0), m_BitmapSize(0) {
		Load(Filename);
	}
	
	~CBitmap() {
		Dispose();
	}
	
	void Dispose() {
		if (m_BitmapData) {
			delete[] m_BitmapData;
			m_BitmapData = 0;
		}
		memset(&m_BitmapFileHeader, 0, sizeof(m_BitmapFileHeader));
		memset(&m_BitmapHeader, 0, sizeof(m_BitmapHeader));
	}
	
	/* Load specified Bitmap and stores it as RGBA in an internal buffer */
	
	bool Load(const char *Filename) {
		std::ifstream file(Filename, std::ios::binary | std::ios::in);
		
		if (file.bad()) {
			return false;
		}

		if (file.is_open() == false) {
			return false;
		}
		
		Dispose();
		
		file.read((char*) &m_BitmapFileHeader, BITMAP_FILEHEADER_SIZE);
		if (m_BitmapFileHeader.Signature != BITMAP_SIGNATURE) {
			return false;
		}

		file.read((char*) &m_BitmapHeader, sizeof(BITMAP_HEADER_SHORT));
		
		/* Load Color Table */
		
		file.seekg(BITMAP_FILEHEADER_SIZE + m_BitmapHeader.HeaderSize, std::ios::beg);

		unsigned int ColorTableSize = 256;
		
		// Always allocate full sized color table

		uint8_t* ColorTable = new uint8_t[ColorTableSize]; // std::bad_alloc exception should be thrown if memory is not available
		file.read((char*) ColorTable, m_BitmapHeader.ClrUsed);

		/* ... Color Table for 16 bits images are not supported yet */	
		
		m_BitmapSize = GetWidth() * GetHeight();
		m_BitmapData = new uint8_t[m_BitmapSize];
		
		unsigned int LineWidth = ((GetWidth() * GetBitCount() / 8) + 3) & ~3;
		uint8_t *Line = new uint8_t[LineWidth];
		
		file.seekg(m_BitmapFileHeader.BitsOffset, std::ios::beg);

		//std::cout << LineWidth << " " <<  m_BitmapHeader.ClrUsed << " " << m_BitmapFileHeader.BitsOffset << std::endl;

		// Ugly fix to BMP being stored with bottom left first.
		int Index = GetWidth() * GetHeight();
		//int Index = 0;
		bool Result = true;

		for (unsigned int i = 0; i < GetHeight(); i++) {
			file.read((char*) Line, LineWidth);
			Index -= GetWidth();

			uint8_t *LinePtr = Line;
				
			for (unsigned int j = 0; j < GetWidth(); j++) {
				m_BitmapData[Index] = *((uint8_t*) LinePtr);
				Index++;
				LinePtr++;
			}

			Index -= GetWidth();
		}
		
		delete [] Line;

		file.close();
		return Result;
	}

	bool Save(const char* Filename) {
		bool Result = true;

		std::ofstream file(Filename, std::ios::out | std::ios::binary);

		if (file.is_open() == false) {
			return false;
		}

		BITMAP_FILEHEADER bfh;
		BITMAP_HEADER_SHORT bh;
		memset(&bfh, 0, sizeof(bfh));
		memset(&bh, 0, sizeof(bh));

		bfh.Signature = BITMAP_SIGNATURE;
		bfh.BitsOffset = BITMAP_FILEHEADER_SIZE + sizeof(BITMAP_HEADER_SHORT) + 256 * 4;
		bfh.Size = GetWidth() * GetHeight() + bfh.BitsOffset;
	
		bh.HeaderSize = sizeof(BITMAP_HEADER_SHORT);
		bh.BitCount = 8;
		bh.Compression = 0;  // No compression for grayscale.

		unsigned int LineWidth = (GetWidth() + 3) & ~3;

		bh.Planes = 1;
		bh.Height = GetHeight();
		bh.Width = GetWidth();
		bh.SizeImage = (LineWidth * 8 * GetHeight()) / 8;
		bh.PelsPerMeterX = 3780;
		bh.PelsPerMeterY = 3780;
		
		// Create an 8-bit Grayscale palette; i.e. all the same.
		uint8_t* Palette = new uint8_t[256 * 4];
		for (auto i=0; i < 256 * 4; i++)
			Palette[i] = i / 4;

		uint8_t* Bitmap = new uint8_t[bh.SizeImage];
		GetBits(Bitmap);
	
		file.write((char*) &bfh, sizeof(BITMAP_FILEHEADER));
		file.write((char*) &bh, sizeof(BITMAP_HEADER_SHORT));
		file.write((char*) Palette, 256 * 4);

		// BMPs are stored upside down for whatever reason...
		for (int row=GetHeight()-1; row>=0; row--) {
			file.write((char*) Bitmap + row * LineWidth, LineWidth);
		}
		
				
		delete [] Bitmap;

		file.close();
		return Result;
	}

	unsigned int GetWidth() {
		/* Add plausibility test */
		// if (abs(m_BitmapHeader.Width) > 8192) {
		//	m_BitmapHeader.Width = 8192;
		// }
		return m_BitmapHeader.Width < 0 ? -m_BitmapHeader.Width : m_BitmapHeader.Width;
	}
	
	unsigned int GetHeight() {
		/* Add plausibility test */
		// if (abs(m_BitmapHeader.Height) > 8192) {
		//	m_BitmapHeader.Height = 8192;
		// }
		return m_BitmapHeader.Height < 0 ? -m_BitmapHeader.Height : m_BitmapHeader.Height;
	}
	
	unsigned int GetBitCount() {
		/* Add plausibility test */
		// if (m_BitmapHeader.BitCount > 32) {
		//	m_BitmapHeader.BitCount = 32;
		// }
		return m_BitmapHeader.BitCount;
	}
	
	/* Copies internal RGBA buffer to user specified buffer */
	void GetBits(uint8_t* Buffer) {
		unsigned int Size = m_BitmapSize;
		memcpy(Buffer, m_BitmapData, Size);
	}

	/* Set Bitmap Bits from scratch. Removes any data that is in the image. */
	void SetBits(uint8_t* Buffer, unsigned int Width, unsigned int Height) {
		if (Buffer == 0) {
			return;
		}

		Dispose();

		uint8_t *BufferPtr = (uint8_t*) Buffer;
		m_BitmapHeader.Width = Width;
		m_BitmapHeader.Height = Height;
		m_BitmapHeader.BitCount = 8;
		m_BitmapHeader.Compression = 0; 

		m_BitmapSize = GetWidth() * GetHeight();
		m_BitmapData = new uint8_t[m_BitmapSize];

		for (unsigned int i = 0; i < m_BitmapSize; i++) {
			m_BitmapData[i] = *((uint8_t*) BufferPtr);
			BufferPtr++;
		}
	}
};

#endif
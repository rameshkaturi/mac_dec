#include "All.h"
#include "Prepare.h"

namespace APE
{

static uint32 CRC32_TABLE[8][256];

static uint32 CRC_reflect(uint32 ref, char ch)
{
    uint32 value = 0;

    for (int i = 1; i < (ch + 1); ++i)
    {
        if (ref & 1) value |= 1 << (ch - i);
        ref >>= 1;
    }

    return value;
}

static bool CRC_init()
{
    uint32 polynomial = 0x4c11db7;

    for (int i = 0; i <= 0xFF; i++)
    {
        uint32 crc = CRC_reflect(i, 8) << 24;

        for (int j = 0; j < 8; j++)
            crc = (crc << 1) ^ (crc & (1 << 31) ? polynomial : 0);

        CRC32_TABLE[0][i] = CRC_reflect(crc, 32);
    }

    for (int i = 0; i <= 0xFF; i++)
        for (int j = 1; j < 8; j++)
            CRC32_TABLE[j][i] = CRC32_TABLE[0][CRC32_TABLE[j - 1][i] & 0xFF] ^ (CRC32_TABLE[j - 1][i] >> 8);

    return true;
}

static bool CRC_initialized = CRC_init();

static uint32 CRC_update(uint32 crc, const unsigned char * pData, int nBytes)
{
	while (nBytes >= 8)
	{
		crc  ^= pData[3] << 24 | pData[2] << 16 | pData[1] << 8 | pData[0];

		crc   =	CRC32_TABLE[7][ crc	   & 0xFF] ^ CRC32_TABLE[6][(crc >>  8) & 0xFF] ^
			CRC32_TABLE[5][(crc >> 16) & 0xFF] ^ CRC32_TABLE[4][ crc >> 24	      ] ^
			CRC32_TABLE[3][pData[4]		 ] ^ CRC32_TABLE[2][pData[5]	      ] ^
			CRC32_TABLE[1][pData[6]		 ] ^ CRC32_TABLE[0][pData[7]	      ];

		pData += 8;
		nBytes -= 8;
	}

	while (nBytes--) crc = (crc >> 8) ^ CRC32_TABLE[0][(crc & 0xFF) ^ *pData++];

	return crc;
}

int CPrepare::Prepare(const unsigned char * pRawData, int nBytes, const WAVEFORMATEX * pWaveFormatEx, int * pOutputX, int *pOutputY, unsigned int *pCRC, int *pSpecialCodes, intn *pPeakLevel)
{
    // error check the parameters
    if (pRawData == NULL || pWaveFormatEx == NULL)
        return ERROR_BAD_PARAMETER;

    // initialize the pointers that got passed in
    *pCRC = 0xFFFFFFFF;
    *pSpecialCodes = 0;

    // variables
    uint32 CRC = 0xFFFFFFFF;
    const int nTotalBlocks = nBytes / pWaveFormatEx->nBlockAlign;
    int R,L;

    // calculate CRC
    CRC = CRC_update(CRC, pRawData, nTotalBlocks * pWaveFormatEx->nChannels * (pWaveFormatEx->wBitsPerSample / 8));

    // the prepare code

    if (pWaveFormatEx->wBitsPerSample == 8) 
    {
        if (pWaveFormatEx->nChannels == 2) 
        {
            for (int nBlockIndex = 0; nBlockIndex < nTotalBlocks; nBlockIndex++) 
            {
                R = (int) (*((unsigned char *) pRawData++) - 128);
                L = (int) (*((unsigned char *) pRawData++) - 128);
                
                // check the peak
                if (labs(L) > *pPeakLevel)
                    *pPeakLevel = labs(L);
                if (labs(R) > *pPeakLevel)
                    *pPeakLevel = labs(R);

                // convert to x,y
                pOutputY[nBlockIndex] = L - R;
                pOutputX[nBlockIndex] = R + (pOutputY[nBlockIndex] / 2);
            }
        }
        else if (pWaveFormatEx->nChannels == 1) 
        {
            for (int nBlockIndex = 0; nBlockIndex < nTotalBlocks; nBlockIndex++) 
            {
                R = (int) (*((unsigned char *) pRawData++) - 128);
                
                // check the peak
                if (labs(R) > *pPeakLevel)
                    *pPeakLevel = labs(R);

                // convert to x,y
                pOutputX[nBlockIndex] = R;
            }
        }
    }
    else if (pWaveFormatEx->wBitsPerSample == 24) 
    {
        if (pWaveFormatEx->nChannels == 2) 
        {
            for (int nBlockIndex = 0; nBlockIndex < nTotalBlocks; nBlockIndex++) 
            {
                uint32 nTemp = 0;
                
                nTemp |= (*pRawData++ << 0);
                nTemp |= (*pRawData++ << 8);
                nTemp |= (*pRawData++ << 16);

                if (nTemp & 0x800000)
                    R = (int) (nTemp & 0x7FFFFF) - 0x800000;
                else
                    R = (int) (nTemp & 0x7FFFFF);

                nTemp = 0;

                nTemp |= (*pRawData++ << 0);                
                nTemp |= (*pRawData++ << 8);                
                nTemp |= (*pRawData++ << 16);
                                
                if (nTemp & 0x800000)
                    L = (int) (nTemp & 0x7FFFFF) - 0x800000;
                else
                    L = (int) (nTemp & 0x7FFFFF);

                // check the peak
                if (labs(L) > *pPeakLevel)
                    *pPeakLevel = labs(L);
                if (labs(R) > *pPeakLevel)
                    *pPeakLevel = labs(R);

                // convert to x,y
                pOutputY[nBlockIndex] = L - R;
                pOutputX[nBlockIndex] = R + (pOutputY[nBlockIndex] / 2);

            }
        }
        else if (pWaveFormatEx->nChannels == 1) 
        {
            for (int nBlockIndex = 0; nBlockIndex < nTotalBlocks; nBlockIndex++) 
            {
                uint32 nTemp = 0;
                
                nTemp |= (*pRawData++ << 0);                
                nTemp |= (*pRawData++ << 8);                
                nTemp |= (*pRawData++ << 16);
                
                if (nTemp & 0x800000)
                    R = (int) (nTemp & 0x7FFFFF) - 0x800000;
                else
                    R = (int) (nTemp & 0x7FFFFF);
    
                // check the peak
                if (labs(R) > *pPeakLevel)
                    *pPeakLevel = labs(R);

                // convert to x,y
                pOutputX[nBlockIndex] = R;
            }
        }
    }
    else 
    {
        if (pWaveFormatEx->nChannels == 2) 
        {
            int LPeak = 0;
            int RPeak = 0;
            int nBlockIndex = 0;
            for (nBlockIndex = 0; nBlockIndex < nTotalBlocks; nBlockIndex++) 
            {
                R = (int) *((int16 *) pRawData); pRawData += 2;
                L = (int) *((int16 *) pRawData); pRawData += 2;

                // check the peak
                if (labs(L) > LPeak)
                    LPeak = labs(L);
                if (labs(R) > RPeak)
                    RPeak = labs(R);

                // convert to x,y
                pOutputY[nBlockIndex] = L - R;
                pOutputX[nBlockIndex] = R + (pOutputY[nBlockIndex] / 2);
            }

            if (LPeak == 0) { *pSpecialCodes |= SPECIAL_FRAME_LEFT_SILENCE; }
            if (RPeak == 0) { *pSpecialCodes |= SPECIAL_FRAME_RIGHT_SILENCE; }
            if (ape_max(LPeak, RPeak) > *pPeakLevel) 
            {
                *pPeakLevel = ape_max(LPeak, RPeak);
            }

            // check for pseudo-stereo files
            nBlockIndex = 0;
            while (pOutputY[nBlockIndex++] == 0) 
            {
                if (nBlockIndex == (nBytes / 4)) 
                {
                    *pSpecialCodes |= SPECIAL_FRAME_PSEUDO_STEREO;
                    break;
                }
            }
        }
        else if (pWaveFormatEx->nChannels == 1) 
        {
            int nPeak = 0;
            for (int nBlockIndex = 0; nBlockIndex < nTotalBlocks; nBlockIndex++) 
            {
                R = (int) *((int16 *) pRawData); pRawData += 2;
                
                // check the peak
                if (labs(R) > nPeak)
                    nPeak = labs(R);

                //convert to x,y
                pOutputX[nBlockIndex] = R;
            }

            if (nPeak > *pPeakLevel)
                *pPeakLevel = nPeak;
            if (nPeak == 0) { *pSpecialCodes |= SPECIAL_FRAME_MONO_SILENCE; }
        }
    }

    CRC = CRC ^ 0xFFFFFFFF;

    // add the special code
    CRC >>= 1;

    if (*pSpecialCodes != 0) 
    {
        CRC |= (1 << 31);
    }

    *pCRC = CRC;

    return ERROR_SUCCESS;
}

void CPrepare::Unprepare(int X, int Y, const WAVEFORMATEX * pWaveFormatEx, unsigned char * pOutput, unsigned int * pCRC)
{
    #define CALCULATE_CRC_BYTE    *pCRC = (*pCRC >> 8) ^ CRC32_TABLE[0][(*pCRC & 0xFF) ^ *pOutput++];
    
    // decompress and convert from (x,y) -> (l,r)
    if (pWaveFormatEx->nChannels == 2) 
    {
        if (pWaveFormatEx->wBitsPerSample == 16) 
        {
            // get the right and left values
            int nR = X - (Y / 2);
            int nL = nR + Y;

            // error check (for overflows)
            if ((nR < -32768) || (nR > 32767) || (nL < -32768) || (nL > 32767))
            {
                throw(-1);
            }

            *(int16 *) pOutput = (int16) nR;
            CALCULATE_CRC_BYTE
            CALCULATE_CRC_BYTE
                
            *(int16 *) pOutput = (int16) nL;
            CALCULATE_CRC_BYTE
            CALCULATE_CRC_BYTE
        }
        else if (pWaveFormatEx->wBitsPerSample == 8) 
        {
            unsigned char R = (X - (Y / 2) + 128);
            *pOutput = R;
            CALCULATE_CRC_BYTE
            *pOutput = (unsigned char) (R + Y);
            CALCULATE_CRC_BYTE
        }
        else if (pWaveFormatEx->wBitsPerSample == 24) 
        {
            int32 RV, LV;

            RV = X - (Y / 2);
            LV = RV + Y;
            
            uint32 nTemp = 0;
            if (RV < 0)
                nTemp = ((uint32) (RV + 0x800000)) | 0x800000;
            else
                nTemp = (uint32) RV;    
            
            *pOutput = (unsigned char) ((nTemp >> 0) & 0xFF);
            CALCULATE_CRC_BYTE
            *pOutput = (unsigned char) ((nTemp >> 8) & 0xFF);
            CALCULATE_CRC_BYTE
            *pOutput = (unsigned char) ((nTemp >> 16) & 0xFF);
            CALCULATE_CRC_BYTE

            nTemp = 0;
            if (LV < 0)
                nTemp = ((uint32) (LV + 0x800000)) | 0x800000;
            else
                nTemp = (uint32) LV;    
            
            *pOutput = (unsigned char) ((nTemp >> 0) & 0xFF);
            CALCULATE_CRC_BYTE
            
            *pOutput = (unsigned char) ((nTemp >> 8) & 0xFF);
            CALCULATE_CRC_BYTE
            
            *pOutput = (unsigned char) ((nTemp >> 16) & 0xFF);
            CALCULATE_CRC_BYTE
        }
    }
    else if (pWaveFormatEx->nChannels == 1) 
    {
        if (pWaveFormatEx->wBitsPerSample == 16) 
        {
            int16 R = X;
                
            *(int16 *) pOutput = (int16) R;
            CALCULATE_CRC_BYTE
            CALCULATE_CRC_BYTE
        }
        else if (pWaveFormatEx->wBitsPerSample == 8) 
        {
            unsigned char R = X + 128;
            *pOutput = R;
            CALCULATE_CRC_BYTE
        }
        else if (pWaveFormatEx->wBitsPerSample == 24) 
        {
            int32 RV = X;
            
            uint32 nTemp = 0;
            if (RV < 0)
                nTemp = ((uint32) (RV + 0x800000)) | 0x800000;
            else
                nTemp = (uint32) RV;    
            
            *pOutput = (unsigned char) ((nTemp >> 0) & 0xFF);
            CALCULATE_CRC_BYTE
            *pOutput = (unsigned char) ((nTemp >> 8) & 0xFF);
            CALCULATE_CRC_BYTE
            *pOutput = (unsigned char) ((nTemp >> 16) & 0xFF);
            CALCULATE_CRC_BYTE
        }
    }
}

#ifdef APE_BACKWARDS_COMPATIBILITY

int CPrepare::UnprepareOld(int *pInputX, int *pInputY, intn nBlocks, const WAVEFORMATEX *pWaveFormatEx, unsigned char *pRawData, unsigned int *pCRC, int *pSpecialCodes, intn nFileVersion)
{
	// the CRC that will be figured during decompression
	uint32 CRC = 0xFFFFFFFF;

	// decompress and convert from (x,y) -> (l,r)
	if (pWaveFormatEx->nChannels == 2) 
	{
		// convert the x,y data to raw data
		if (pWaveFormatEx->wBitsPerSample == 16) 
		{
			int16 R;
			unsigned char *Buffer = &pRawData[0];
			int *pX = pInputX;
			int *pY = pInputY;

			for (; pX < &pInputX[nBlocks]; pX++, pY++) 
			{
				R = *pX - (*pY / 2);

				*(int16 *) Buffer = (int16) R;
				CRC = (CRC >> 8) ^ CRC32_TABLE[0][(CRC & 0xFF) ^ *Buffer++];
				CRC = (CRC >> 8) ^ CRC32_TABLE[0][(CRC & 0xFF) ^ *Buffer++];

				*(int16 *) Buffer = (int16) R + *pY;
				CRC = (CRC >> 8) ^ CRC32_TABLE[0][(CRC & 0xFF) ^ *Buffer++];
				CRC = (CRC >> 8) ^ CRC32_TABLE[0][(CRC & 0xFF) ^ *Buffer++];
			}
		}
		else if (pWaveFormatEx->wBitsPerSample == 8) 
		{
			unsigned char *R = (unsigned char *) &pRawData[0];
			unsigned char *L = (unsigned char *) &pRawData[1];

			if (nFileVersion > 3830) 
			{
				for (int SampleIndex = 0; SampleIndex < nBlocks; SampleIndex++, L+=2, R+=2) 
				{
					*R = (unsigned char) (pInputX[SampleIndex] - (pInputY[SampleIndex] / 2) + 128);
					CRC = (CRC >> 8) ^ CRC32_TABLE[0][(CRC & 0xFF) ^ *R];
					*L = (unsigned char) (*R + pInputY[SampleIndex]);
					CRC = (CRC >> 8) ^ CRC32_TABLE[0][(CRC & 0xFF) ^ *L];
				}
			}
			else 
			{
				for (int SampleIndex = 0; SampleIndex < nBlocks; SampleIndex++, L+=2, R+=2)
				{
					*R = (unsigned char) (pInputX[SampleIndex] - (pInputY[SampleIndex] / 2));
					CRC = (CRC >> 8) ^ CRC32_TABLE[0][(CRC & 0xFF) ^ *R];
					*L = (unsigned char) (*R + pInputY[SampleIndex]);
					CRC = (CRC >> 8) ^ CRC32_TABLE[0][(CRC & 0xFF) ^ *L];
				}
			}
		}
		else if (pWaveFormatEx->wBitsPerSample == 24) 
		{
			unsigned char *Buffer = (unsigned char *) &pRawData[0];
			int32 RV, LV;

			for (int SampleIndex = 0; SampleIndex < nBlocks; SampleIndex++)
			{
				RV = pInputX[SampleIndex] - (pInputY[SampleIndex] / 2);
				LV = RV + pInputY[SampleIndex];

				uint32 nTemp = 0;
				if (RV < 0)
					nTemp = ((uint32) (RV + 0x800000)) | 0x800000;
				else
					nTemp = (uint32) RV;    

				*Buffer = (unsigned char) ((nTemp >> 0) & 0xFF);
				CRC = (CRC >> 8) ^ CRC32_TABLE[0][(CRC & 0xFF) ^ *Buffer++];

				*Buffer = (unsigned char) ((nTemp >> 8) & 0xFF);
				CRC = (CRC >> 8) ^ CRC32_TABLE[0][(CRC & 0xFF) ^ *Buffer++];

				*Buffer = (unsigned char) ((nTemp >> 16) & 0xFF);
				CRC = (CRC >> 8) ^ CRC32_TABLE[0][(CRC & 0xFF) ^ *Buffer++];

				nTemp = 0;
				if (LV < 0)
					nTemp = ((uint32) (LV + 0x800000)) | 0x800000;
				else
					nTemp = (uint32) LV;    

				*Buffer = (unsigned char) ((nTemp >> 0) & 0xFF);
				CRC = (CRC >> 8) ^ CRC32_TABLE[0][(CRC & 0xFF) ^ *Buffer++];

				*Buffer = (unsigned char) ((nTemp >> 8) & 0xFF);
				CRC = (CRC >> 8) ^ CRC32_TABLE[0][(CRC & 0xFF) ^ *Buffer++];

				*Buffer = (unsigned char) ((nTemp >> 16) & 0xFF);
				CRC = (CRC >> 8) ^ CRC32_TABLE[0][(CRC & 0xFF) ^ *Buffer++];
			}
		}
	}
	else if (pWaveFormatEx->nChannels == 1) 
	{
		// convert to raw data
		if (pWaveFormatEx->wBitsPerSample == 8) 
		{
			unsigned char *R = (unsigned char *) &pRawData[0];

			if (nFileVersion > 3830) 
			{
				for (int SampleIndex = 0; SampleIndex < nBlocks; SampleIndex++, R++)
				{
					*R = pInputX[SampleIndex] + 128;
					CRC = (CRC >> 8) ^ CRC32_TABLE[0][(CRC & 0xFF) ^ *R];
				}
			}
			else 
			{
				for (int SampleIndex = 0; SampleIndex < nBlocks; SampleIndex++, R++)
				{
					*R = (unsigned char) (pInputX[SampleIndex]);
					CRC = (CRC >> 8) ^ CRC32_TABLE[0][(CRC & 0xFF) ^ *R];
				}
			}
		}
		else if (pWaveFormatEx->wBitsPerSample == 24) 
		{
			unsigned char *Buffer = (unsigned char *) &pRawData[0];
			int32 RV;
			for (int SampleIndex = 0; SampleIndex<nBlocks; SampleIndex++) 
			{
				RV = pInputX[SampleIndex];

				uint32 nTemp = 0;
				if (RV < 0)
					nTemp = ((uint32) (RV + 0x800000)) | 0x800000;
				else
					nTemp = (uint32) RV;    

				*Buffer = (unsigned char) ((nTemp >> 0) & 0xFF);
				CRC = (CRC >> 8) ^ CRC32_TABLE[0][(CRC & 0xFF) ^ *Buffer++];

				*Buffer = (unsigned char) ((nTemp >> 8) & 0xFF);
				CRC = (CRC >> 8) ^ CRC32_TABLE[0][(CRC & 0xFF) ^ *Buffer++];

				*Buffer = (unsigned char) ((nTemp >> 16) & 0xFF);
				CRC = (CRC >> 8) ^ CRC32_TABLE[0][(CRC & 0xFF) ^ *Buffer++];
			}
		}
		else 
		{
			unsigned char *Buffer = &pRawData[0];

			for (int SampleIndex = 0; SampleIndex < nBlocks; SampleIndex++) 
			{
				*(int16 *) Buffer = (int16) (pInputX[SampleIndex]);
				CRC = (CRC >> 8) ^ CRC32_TABLE[0][(CRC & 0xFF) ^ *Buffer++];
				CRC = (CRC >> 8) ^ CRC32_TABLE[0][(CRC & 0xFF) ^ *Buffer++];
			}
		}
	}

	CRC = CRC ^ 0xFFFFFFFF;

	*pCRC = CRC;

	return 0;
}

#endif // #ifdef APE_BACKWARDS_COMPATIBILITY

}
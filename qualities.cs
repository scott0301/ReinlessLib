using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ReinlessLib
{
    public static partial class CReinlessLib
    {
        //***********************************************************************************
        // Image Quality measurement 
        //***********************************************************************************

        #region QUALITY MEASUREMENT
        public static double /***/HC_QUALITY_GetQuality_MSE(byte[] rawSrc, int imageW1, int imageH1, byte[] rawDest, int imageW2, int imageH2)
        {
            if (imageW1 != imageW2 || imageH1 != imageH2)
                return -1;

            byte[] rawImage1 = rawSrc;
            byte[] rawImage2 = rawDest;

            double fDiffSum = 0.0;

            for (int y = 0; y < imageH1; y++)
            {
                for (int x = 0; x < imageW1; x++)
                {
                    fDiffSum += Math.Pow(Convert.ToDouble(rawImage1[y * imageW1 + x]) - Convert.ToDouble(rawImage2[y * imageW2 + x]), 2);
                }
            }

            return fDiffSum / Convert.ToDouble(imageW1 * imageH1);
        }
        public static double /***/HC_QUALITY_GetQuality_NAE(byte[] rawSrc, int imageW1, int imageH1, byte[] rawDest, int imageW2, int imageH2)
        {
            if (imageW1 != imageW2 || imageH1 != imageH2)
                return -1;

            byte[] rawImage1 = rawSrc;
            byte[] rawImage2 = rawDest;

            double fSumDiff = 0.0;
            double fSumOrigin = 0.0;

            for (int y = 0; y < imageH1; y++)
            {
                for (int x = 0; x < imageW1; x++)
                {
                    fSumOrigin += rawImage1[y * imageW1 + x];
                    fSumDiff += Math.Abs(Convert.ToDouble(rawImage1[y * imageW1 + x]) - Convert.ToDouble(rawImage2[y * imageW2 + x]));
                }
            }

            return fSumDiff / fSumOrigin;
        }
        public static double /***/HC_QUALITY_GetQuality_NCC(byte[] rawSrc, int imageW1, int imageH1, byte[] rawDest, int imageW2, int imageH2)
        {
            if (imageW1 != imageW2 || imageH1 != imageH2)
                return -1;

            byte[] rawImage1 = rawSrc;
            byte[] rawImage2 = rawDest;

            double fSquareOrigin = 0.0;
            double fSquareEach = 0.0;

            for (int y = 0; y < imageH1; y++)
            {
                for (int x = 0; x < imageW1; x++)
                {
                    fSquareEach += rawImage1[y * imageW1 + x] * rawImage2[y * imageW2 + x];
                    fSquareOrigin += rawImage1[y * imageW1 + x] * rawImage1[y * imageW1 + x];
                }
            }

            return fSquareEach / fSquareOrigin;
        }
        public static double /***/HC_QUALITY_GetQuality_PSNR(byte[] rawSrc, int imageW1, int imageH1, byte[] rawDest, int imageW2, int imageH2)
        {
            double MSE = HC_QUALITY_GetQuality_MSE(rawSrc, imageW1, imageH1, rawDest, imageW2, imageH2);
            double PSNR = 99.0;

            if (MSE > 0)
            {
                PSNR = 10 * Math.Log(255.0 * 255.0 / MSE) / Math.Log(10);
            }
            return PSNR;
        }
        public static void /*****/HC_QUALITY_GetQuality_CMN(byte[] rawSrc, int imageW1, int imageH1, byte[] rawDest, int imageW2, int imageH2, out double fMin, out double fMax, out double fAVG)
        {

            fMin = 999.0;
            fMax = 0x00;
            fAVG = 0x00;

            if (imageW1 != imageW2 || imageH1 != imageH2)
                return;

            byte[] rawImage1 = rawSrc;
            byte[] rawImage2 = rawDest;

            double fDiff = 0.0;

            for (int y = 0; y < imageH1; y++)
            {
                for (int x = 0; x < imageW1; x++)
                {
                    fDiff = rawImage1[y * imageW1 + x] - rawImage2[y * imageW2 + x];

                    if (fMax < Math.Abs(fDiff)) fMax = Math.Abs(fDiff);
                    if (fMin > Math.Abs(fDiff)) fMin = Math.Abs(fDiff);

                    fAVG += fDiff;
                }
            }
            fAVG /= imageW1 * imageH1;

        }
        public static double /***/HC_QUALITY_GetQuality_SelfNoise(byte[] rawImage, int imageW, int imageH)
        {
            double fSigma = 0;

            double fBaseValue = 0;

            for (int y = 1; y < imageH - 1; y++)
            {
                for (int x = 1; x < imageW - 1; x++)
                {
                    fBaseValue += Math.Abs((rawImage[y * imageW + x - 1] * (+1)) + (rawImage[(y - 1) * imageW + x] * (-2)) + (rawImage[y * imageW + x + 1] * (+1)) +
                                           (rawImage[y * imageW + x - 1] * (-2)) + (rawImage[(y + 0) * imageW + x] * (+4)) + (rawImage[y * imageW + x + 1] * (-2)) +
                                           (rawImage[y * imageW + x - 1] * (+1)) + (rawImage[(y + 1) * imageW + x] * (-2)) + (rawImage[y * imageW + x + 1] * (+1)));
                }
            }

            fSigma = fBaseValue * Math.Sqrt(0.5 * Math.PI) / (6 * (imageW - 2) * (imageH - 2));

            return fSigma;
        }
        #endregion

       
    }
}

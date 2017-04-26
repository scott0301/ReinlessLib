using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Drawing;
using System.Drawing.Imaging;

namespace ReinlessLib
{
    public static partial class CReinlessLib
    {
        public static byte[] HC_ARITH_ADD(byte[] rawImage1, byte[] rawImage2, int imageW, int imageH)
        {
            byte[] newImage = new byte[imageW * imageH];

            for (int i = 0; i < imageW * imageH; i++)
            {
                int nValue = rawImage1[i] + rawImage2[i];
                newImage[i] = nValue > 255 ? (byte)255 : (byte)nValue;
            }

            return newImage;
        }
        public static byte[] HC_ARITH_MID(byte[] rawImage1, byte[] rawImage2, int imageW, int imageH)
        {
            byte[] newImage = new byte[imageW * imageH];

            for (int i = 0; i < imageW * imageH; i++)
            {
                int nValue = (rawImage1[i] + rawImage2[i]) /2;
                newImage[i] = nValue > 255 ? (byte)255 : (byte)nValue;
            }

            return newImage;
        }
        public static byte[] HC_ARITH_SUB(byte[] rawImage1, byte[] rawImage2, int imageW, int imageH)
        {
            byte[] newImage = new byte[imageW * imageH];

            for (int i = 0; i < imageW * imageH; i++)
            {
                int nValue = rawImage1[i] - rawImage2[i];
                newImage[i] = nValue < 0 ? (byte)0 : (byte)nValue;
            }
            return newImage;
        }
        public static byte[] HC_ARITH_MUL(byte[] rawImage, int imageW, int imageH, double fValue)
        {
            double [] fImage = HC_CONV_Byte2Double( rawImage);

            Parallel.For(0, imageW * imageH, i =>
            {
                double fPixel = rawImage[i] * fValue;
                rawImage[i] = fPixel > 255 ? (byte)255 : (byte)fPixel;
            });

            //HC_HISTO_AvoidSaturation(ref fImage);
            //
            //rawImage = HC_CONV_Double2Byte(fImage);

            return rawImage;
        }
        public static byte[] HC_ARITH_AND(byte[] rawImage1, byte []rawImage2, int imageW, int imageH)
        {
            byte[] newImage = new byte[imageW * imageH];

            for (int i = 0; i < imageW * imageH; i++)
            {
                if (rawImage1[i] != 0 && rawImage2[i] != 0)
                {
                    newImage[i] = (byte)255;
                }
            }
            return newImage;
        }
        public static byte[] HC_ARITH_OR(byte [] rawImage1, byte[] rawImage2, int imageW, int imageH)
        {
            byte [] newImage = new byte[imageW*imageH];

            for(int i = 0; i < imageW*imageH; i++)
            {
                if( rawImage1[i] != 0 || rawImage2[i] != 0)
                {
                    newImage[i] = (byte)255;
                }
            }
            return newImage;
        }
        public static byte[] HC_ARITH_XOR(byte[] rawImage1, byte [] rawImage2, int imageW, int imageH)
        {
            byte[] newImage = new byte[imageW * imageH];

            for (int i = 0; i < imageW * imageH; i++)
            {
                int nValue1 = rawImage1[i];
                int nValue2 = rawImage1[i];

                newImage[i] = rawImage1[i] == rawImage2[i] ? rawImage1[i] : (byte)255;
            }
            return newImage;
        }

        public static byte[] /**/HC_CONV_BlendedImage(byte[] i1, byte[] i2, int imageW, int imageH, int nBlend)
        {
            byte[] returnRaw = new byte[imageW * imageH];

            //for (int i = 0; i < returnRaw.Length; i++)
            Parallel.For(0, returnRaw.Length, i =>
            {
                double fValue = (i1[i] * ((100 - nBlend) / 100.0)) + (i2[i] * (nBlend / 100.0));
                returnRaw[i] = fValue >= 255 ? (byte)255 : fValue < 0 ? (byte)0 : Convert.ToByte(fValue);

            });

            return returnRaw;
        }

        
    }
}

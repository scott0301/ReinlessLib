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
        //******************************************************************************************
        // Image Conversion
        //******************************************************************************************

        #region IMAGE CONVERSION
        public static void /*******/HC_CONV_Double2Byte(double[] fArray, byte[] byteArray)
        {
            if (fArray.Length == byteArray.Length)
            {
                Parallel.For(0, fArray.Length, i =>
                {
                    double fValue = fArray[i];
                    fValue = fValue < 0x0 ? 0x0 : fValue > 0xff ? 0xff : fValue;
                    byteArray[i] = (byte)fValue;
                });

            }
        }
        public static byte[] /*****/HC_CONV_Double2Byte(double[] fArray)
        {
            byte[] rawImage = new byte[fArray.Length];

            Parallel.For(0, rawImage.Length, i =>
            {
                double fValue = fArray[i];
                fValue = fValue < 0x0 ? 0x0 : fValue > 0xff ? 0xff : fValue;
                rawImage[i] = (byte)fValue;

            });

            return rawImage;
        }
        public static double[] /***/HC_CONV_Byte2Double(byte[] byteArray)
        {
            double[] fArray = new double[byteArray.Length];

            if (fArray.Length == byteArray.Length)
            {
                Parallel.For(0, fArray.Length, i =>
                {
                    fArray[i] = byteArray[i];
                });
            }

            return fArray;
        }
        public static Bitmap/******/HC_CONV_Byte2Bmp(byte[] rawImage, int imageW, int imageH)
        {
            if (imageW == 0 || imageH == 0)
            {
                return new Bitmap(444, 444, PixelFormat.Format24bppRgb);
            }

            Bitmap bmpImage = new Bitmap(imageW, imageH, PixelFormat.Format24bppRgb);

            int nStride = 0, bmpLength = 0;
            byte[] rawBmp = null;

            BitmapData bitmapData = bmpImage.LockBits(new Rectangle(0, 0, imageW, imageH), ImageLockMode.ReadOnly, PixelFormat.Format24bppRgb);
            {
                nStride = Math.Abs(bitmapData.Stride);
                bmpLength = nStride * imageH;

            }
            bmpImage.UnlockBits(bitmapData);


            rawBmp = new byte[bmpLength];

            Parallel.For(0, imageH, y =>
            {
                for (int x = 0; x < imageW; x++)
                {
                    //rawBmp[(y * nStride) + x ] = rawImage[y * imageW + x];
                    rawBmp[(y * nStride) + (x * 3) + 0] =
                    rawBmp[(y * nStride) + (x * 3) + 1] =
                    rawBmp[(y * nStride) + (x * 3) + 2] = rawImage[y * imageW + x];
                }
            });


            bitmapData = bmpImage.LockBits(new Rectangle(0, 0, imageW, imageH), ImageLockMode.ReadWrite, PixelFormat.Format24bppRgb);
            {
                System.Runtime.InteropServices.Marshal.Copy(rawBmp, 0, bitmapData.Scan0, bmpLength);
            }
            bmpImage.UnlockBits(bitmapData);

            return bmpImage;
        }
        public static byte[] /*****/HC_CONV_Bmp2Byte(System.Drawing.Bitmap bmpImage, ref int imageW, ref int imageH)
        {
            imageW = bmpImage.Width;
            imageH = bmpImage.Height;

            int nRealW = 0, nStride = 0, bmpLength = 0;
            byte[] rawBmp = null;

            BitmapData bitmapData = bmpImage.LockBits(new Rectangle(0, 0, imageW, imageH), System.Drawing.Imaging.ImageLockMode.ReadOnly, PixelFormat.Format24bppRgb);
            {
                imageW = bitmapData.Width;
                imageH = bitmapData.Height;
                nRealW = imageW;

                nStride = Math.Abs(bitmapData.Stride);
                bmpLength = nStride * imageH;

                
                rawBmp = new byte[bmpLength];
                System.Runtime.InteropServices.Marshal.Copy(bitmapData.Scan0, rawBmp, 0, bmpLength);
            }
            bmpImage.UnlockBits(bitmapData);

            int nImageW = imageW;
            int nImageH = imageH;

            byte[] rawImage = new byte[imageW * imageH];

            Parallel.For(0, imageH, y =>
            {
                for (int x = 0; x < nImageW; x++)
                {
                    rawImage[y * nImageW + x] = (byte)((rawBmp[(y * nStride) + (x * 3) + 0] + rawBmp[(y * nStride) + (x * 3) + 1] + rawBmp[(y * nStride) + (x * 3) + 2]) / 3);
                }
            });
            return rawImage;
        }
        public static double[] /***/HC_CONV_Bmp2Double(Bitmap bmp, ref int imageW, ref int imageH)
        {
            imageW = bmp.Width;
            imageH = bmp.Height;

            int nRealW = 0, nStride = 0, bmpLength = 0;
            byte[] rawBmp = null;

            BitmapData bitmapData = bmp.LockBits(new Rectangle(0, 0, imageW, imageH), System.Drawing.Imaging.ImageLockMode.ReadOnly, PixelFormat.Format24bppRgb);
            {
                imageW = bitmapData.Width;
                imageH = bitmapData.Height;
                nRealW = imageW;

                nStride = Math.Abs(bitmapData.Stride);
                bmpLength = nStride * imageH;

                rawBmp = new byte[bmpLength];
                System.Runtime.InteropServices.Marshal.Copy(bitmapData.Scan0, rawBmp, 0, bmpLength);
            }
            bmp.UnlockBits(bitmapData);

            int nImageW = imageW;
            int nImageH = imageH;

            double[] fImage = new double[imageW * imageH];

            Parallel.For(0, imageH, y =>
            {
                for (int x = 0; x < nImageW; x++)
                {
                    fImage[y * nImageW + x] = (double)((rawBmp[(y * nStride) + (x * 3) + 0] + rawBmp[(y * nStride) + (x * 3) + 1] + rawBmp[(y * nStride) + (x * 3) + 2]) / 3.0);
                }
            });
            return fImage;
        }
        #endregion

        public static byte[] /*****/HC_CONV_GetNormalizedImage(double[] fImage)
        {
            double MIN = fImage.Min();
            double MAX = fImage.Max();

            if (double.IsNaN(MIN))
                MIN = 0;

            double RANGE = MAX - MIN;

            byte[] rawImage = new byte[fImage.Length];
            
            Parallel.For(0, fImage.Length, i =>
            {
                double fValue = ((fImage[i] - MIN) / (RANGE)) * 255.0 ;

                fValue = double.IsNaN(fValue) == true ? 0 : Math.Floor(fValue);

                rawImage[i] = (byte)fValue;
            });
            return rawImage;
        }
        public static double[] /***/HC_CONV_GetMeanImage(double[] rawImage, int imageW, int imageH, int KernerSize)
        {
            int nKernel = KernerSize;
            int nKH = nKernel / 2;
            int nGap = nKH;

            double[] pad = (double[])ARRAY_Padding_ALL(rawImage, imageW, imageH, nGap);
            double[] mean = new double[imageW * imageH];

            int nCount = nKernel * nKernel;
            double fBlockSum = 0;

            for (int y = nKH; y < imageH + nKH; y++)
            {
                for (int x = nKH; x < imageW + nKH; x++)
                {
                    fBlockSum = 0x00;

                    for (int yy = y - nKH; yy <= y + nKH; yy++)
                    {
                        for (int xx = x - nKH; xx <= x + nKH; xx++)
                        {
                            fBlockSum += pad[yy * (imageW + nGap * 2) + xx];
                        }
                    }

                    mean[(y - nKH) * imageW + (x - nKH)] = fBlockSum / nCount;
                }
            }

            return mean;
        }
        public static double[] /***/HC_CONV_GetMeanImage(byte[] rawImage, int imageW, int imageH, int KernerSize)
        {
            int nKernel = KernerSize;
            int nKH = nKernel / 2;
            int nGap = nKH;

            byte[] pad = (byte[])ARRAY_Padding_ALL(rawImage, imageW, imageH, nGap);
            double[] mean = new double[imageW * imageH];

            int nCount = nKernel * nKernel;
            double fBlockSum = 0;
            int nCheck = 0;
            for (int y = nKH; y < imageH + nKH; y++)
            {
                for (int x = nKH; x < imageW + nKH; x++)
                {
                    fBlockSum = 0x00;
                    nCheck = 0;
                    for (int yy = y - nKH; yy <= y + nKH; yy++)
                    {
                        for (int xx = x - nKH; xx <= x + nKH; xx++)
                        {
                            fBlockSum += pad[yy * (imageW + nGap * 2) + xx];
                            nCheck++;
                        }
                    }

                    mean[(y - nKH) * imageW + (x - nKH)] = fBlockSum / nCount;
                }
            }

            return mean;
        }
        public static double[] /***/HC_CONV_GetPowImage(byte[] rawImage)
        {
            double[] powImage = new double[rawImage.Length];

            Parallel.For(0, rawImage.Length, i =>
            {
                powImage[i] = rawImage[i] * rawImage[i];
            });
            return powImage;
        }
        public static double[] /***/HC_CONV_GetPowImage(double[] rawImage)
        {
            double[] powImage = new double[rawImage.Length];

            Parallel.For(0, rawImage.Length, i =>
            {
                powImage[i] = rawImage[i] * rawImage[i];
            });
            return powImage;
        }

        public static byte EnsureByte(object value)
        {
            if (double.IsNaN((double)value) == true) return 0;
 
            double result /*****/= (double)value; ;
            return (byte)result > 255 ? (byte)255 : (byte)result < 0 ? (byte)0 : (byte)result;
        }

      

        //******************************************************************************************
        // Color Conversion
        //******************************************************************************************
        
        #region COLOR CONVERSION 
        public static Image HC_CONV_Raw2Color(byte[] rawImage, int imageW, int imageH)
        {
            byte[] rawImageC = new byte[imageW * imageH * 3];

            int nMin = rawImage.Min();
            int nMax = rawImage.Max();

            Bitmap bmpImage = new Bitmap(imageW, imageH, PixelFormat.Format24bppRgb);

            int nStride = 0, bmpLength = 0;
            byte[] rawBmp = null;

            BitmapData bitmapData = bmpImage.LockBits(new Rectangle(0, 0, imageW, imageH), ImageLockMode.ReadOnly, PixelFormat.Format24bppRgb);
            {
                nStride = Math.Abs(bitmapData.Stride);
                bmpLength = nStride * imageH;
            }
            bmpImage.UnlockBits(bitmapData);

            rawBmp = new byte[bmpLength];

            Parallel.For(0, imageH, y =>
            {
                for (int x = 0; x < imageW; x++)
                {
                    byte value = rawImage[y * imageW + x];
                    byte[] c = getPseudoColorRaw(value, nMin, nMax);

                    rawBmp[(y * nStride) + (x * 3) + 0] = c[0];
                    rawBmp[(y * nStride) + (x * 3) + 1] = c[1];
                    rawBmp[(y * nStride) + (x * 3) + 2] = c[2];
                }
            });

            bitmapData = bmpImage.LockBits(new Rectangle(0, 0, imageW, imageH), ImageLockMode.ReadWrite, PixelFormat.Format24bppRgb);
            {
                System.Runtime.InteropServices.Marshal.Copy(rawBmp, 0, bitmapData.Scan0, bmpLength);
            }
            bmpImage.UnlockBits(bitmapData);

            return bmpImage;
        }
        public static byte[] getPseudoColorRaw(double value, double min, double max)
        {
            int r = 0;
            int g = 0;
            int b = 0;

            if (value < min + 0.25 * (max - min))
            {
                r = 0;
                g = (int)(4.0 * 255.0 * (value - min) / (max - min));
                b = 255;
            }
            else if (value < min + 0.5 * (max - min))
            {
                r = 0;
                g = 255;
                b = (int)(255.0 + 4.0 * 255.0 * (min + 0.25 * (max - min) - value) / (max - min));
            }
            else if (value < min + 0.75 * (max - min))
            {
                r = (int)(4.0 * 255.0 * (value - min - 0.5 * (max - min)) / (max - min));
                g = 255;
                b = 0;
            }
            else
            {
                r = 255;
                g = (int)(255.0 + 4.0 * 255.0 * (min + 0.75 * (max - min) - value) / (max - min));
                b = 0;
            }

            byte br = (byte)(r);
            byte bg = (byte)(g);
            byte bb = (byte)(b);

            //byte[] c = { br, bg, bb };
            byte[] c = { bb, bg, br };
            return c;
        }

        public static byte[] HC_CONV_Bmp2Monochrome(Bitmap bmp)
        {
            byte [] temp = new byte[0];
            return temp;

            // in common case most frequently and efficiently used.
            // ★ Gray = (Red + Green + Blue) / 3

            //Humans perceive green more strongly than red, and red more strongly than blue
            // Because humans do not perceive all colors equally, the “average method” of grayscale conversion is inaccurate

            // ★ Gray = (Red * 0.3 + Green * 0.59 + Blue * 0.11)

            //original ITU-R recommendation (BT.709, specifically) which is the historical precedent.  
            //This formula, sometimes called Luma,

            //★ Gray = (Red * 0.2126 + Green * 0.7152 + Blue * 0.0722)

            // (BT.601), which calls for slightly different coefficients:
            //★ Gray = (Red * 0.299 + Green * 0.587 + Blue * 0.114)

            //Basically, this takes a color and converts it to its least-saturated variant.
            //★ Gray = ( Max(Red, Green, Blue) + Min(Red, Green, Blue) ) / 2


            //★Maximum decomposition:
            //★ Gray = Max(Red, Green, Blue)
            
            // Minimum decomposition:
            // ★Gray = Min(Red, Green, Blue)

            // Custom # of gray shades
            //Notes:
            //-NumberOfShades is a value between 2 and 256
            //-technically, any grayscale algorithm could be used to calculate AverageValue; it simply provides
            // an initial gray value estimate
            //-the "+ 0.5" addition is an optional parameter that imitates rounding the value of an integer
            // conversion; YMMV depending on which programming language you use, as some round automatically

            // ConversionFactor = 255 / (NumberOfShades - 1)
            //AverageValue = (Red + Green + Blue) / 3
            //Gray = Integer((AverageValue / ConversionFactor) + 0.5) * ConversionFactor
            

           
        }

        #endregion
 
    }
}

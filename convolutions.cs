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
        /************************************************************************************/
        // Convolution related functions
        /************************************************************************************/

        #region PADDING
        public static object ARRAY_Padding_LT(object rawImage, int imageW, int imageH, int nGap)
        {
            object ARRAY_DOUBLE = new double[1];
            object ARRAY_BYTE = new byte[1];
            object returnData = null;

            int newW = imageW + nGap;
            int newH = imageH + nGap;

             
            if (rawImage.GetType() == ARRAY_DOUBLE.GetType())
            {
                #region for double
                int size = sizeof(double);
                double[] fArray = (double[])rawImage;
                double[] newArray = new double[newW * newH];

                int orgY = 0;
                int CopyLength = imageW * size;
                for (int y = nGap; y < newH; y++)
                Buffer.BlockCopy(fArray, (orgY++ * imageW) * size, newArray, ((y * newW) + nGap) * size, CopyLength);

                if (nGap > imageW) return newArray;

                Parallel.For(nGap, newH, y =>{for (int x = 0; x <= nGap; x++){newArray[y * newW + nGap - x] = newArray[y * newW + nGap + x];}});

                CopyLength = newW * size;
                Parallel.For(0, nGap + 1, y =>{Buffer.BlockCopy(newArray, ((nGap + y) * (newW)) * size, newArray, ((nGap - y) * (newW)) * size, CopyLength);});
                returnData = newArray;
                #endregion
            }
            else if (rawImage.GetType() == ARRAY_BYTE.GetType())
            {
                #region for byte

                byte[] byteArray = (byte[])rawImage;
                byte[] newArray = new byte[newW * newH];

                int orgY = 0;

                for (int y = nGap; y < newH; y++)
                Buffer.BlockCopy(byteArray, orgY++ * imageW, newArray, y * newW + nGap, imageW);

                if (nGap > imageW) return newArray;

                Parallel.For(nGap, newH, y =>{for (int x = 0; x <= nGap; x++){newArray[y * newW + nGap - x] = newArray[y * newW + nGap + x];}});
                Parallel.For(0, nGap + 1, y =>{Buffer.BlockCopy(newArray, (nGap + y) * newW, newArray, (nGap - y) * newW, newW);});

                returnData = newArray;
                #endregion
            }
            
            return returnData;
        }
        public static object ARRAY_Padding_RB(object rawImage, int imageW, int imageH, int nGap)
        {
            object ARRAY_DOUBLE = new double[1];
            object ARRAY_BYTE = new byte[1];
            object returnData = null;

            int newW = imageW + nGap;
            int newH = imageH + nGap;

            if (rawImage.GetType() == ARRAY_DOUBLE.GetType())
            {
                int size = sizeof(double);
                double[] fArray = (double[])rawImage;
                double[] newArray = new double[newW * newH];

                int copyLength = imageW * sizeof(double);
                Parallel.For(0, imageH, y =>
                {
                    Buffer.BlockCopy(fArray, (y * imageW) * size, newArray, (y * newW) * size, copyLength);
                });

                // right direction copy
                Parallel.For(0, nGap, x =>
                {
                    for (int y = 0; y < imageH; y++)
                    {
                        newArray[y * newW + imageW + x] = newArray[y * newW + imageW - 1 - x];
                    }
                });

                double[] rawPadVert = new double[newW];

                // bottom direction copy - reverse
                copyLength = newW * sizeof(double);
                Parallel.For(0, nGap, y =>
                {
                    Buffer.BlockCopy(newArray, ((imageH - 1 - y) * newW) * size, rawPadVert, 0, copyLength);
                    Buffer.BlockCopy(rawPadVert, 0, newArray, ((imageH + y) * newW) * size, copyLength);
                });
                returnData = newArray;
            }
            else if (rawImage.GetType() == ARRAY_BYTE.GetType())
            {
                byte[] byteArray = (byte[])rawImage;
                byte[] newArray = new byte[newW * newH];

                Parallel.For(0, imageH, y =>
                {
                    Buffer.BlockCopy(byteArray, y * imageW, newArray, y * newW, imageW);
                });

                // right direction copy
                Parallel.For(0, nGap, x =>
                {
                    for (int y = 0; y < imageH; y++)
                    {
                        newArray[y * newW + imageW + x] = newArray[y * newW + imageW - 1 - x];
                    }
                });

                byte[] rawPadVert = new byte[newW];

                // bottom direction copy - reverse
                Parallel.For(0, nGap, y =>
                {
                    Buffer.BlockCopy(newArray, (imageH - 1 - y) * newW, rawPadVert, 0, newW);
                    Buffer.BlockCopy(rawPadVert, 0, newArray, (imageH + y) * newW, newW);
                });
                returnData = newArray;
            }
            return returnData;
        }
        public static Object ARRAY_Padding_ALL(Object obArray, int arrW, int arrH, int nGap)
        {
            object firstPadding = null;
            object seconPadding = null;
            Parallel.Invoke(()=>{firstPadding = ARRAY_Padding_LT(obArray, arrW, arrH, nGap); });
            Parallel.Invoke(()=>{seconPadding = ARRAY_Padding_RB(firstPadding, arrW + nGap, arrH + nGap, nGap); });
            return seconPadding;
        }
        #endregion

        //*****************************************************************************************
        // Convolution
        //*****************************************************************************************
        
        #region CONVOLUTION
        public static byte [] /*****/HC_FILTER_Convolution(double[] fKernel, byte[] rawImage, int imageW, int imageH)
        {
            double[] fImage = new double[imageW * imageH];
            
            int KSIZE = (int)Math.Sqrt(fKernel.Length);
            int GAP = KSIZE / 2;

            byte[] rawExpanded = (byte[])ARRAY_Padding_ALL(rawImage, imageW, imageH, GAP);


            int imageNewW = imageW + GAP * 2;
            int imageNewH = imageH + GAP * 2;

            //for (int y = GAP; y < imageNewH - GAP; y++)
            Parallel.For(GAP, imageNewH - GAP, y =>
            {
                for (int x = GAP; x < imageNewW - GAP; x++)
                {
                    double kernelSum = 0;
                    for (int j = -GAP; j <= GAP; j++)
                    {
                        for (int k = -GAP; k <= GAP; k++)
                        {
                            kernelSum += (fKernel[(j + GAP) * KSIZE + k + GAP] * rawExpanded[(y - j) * imageNewW + (x - k)]);
                        }
                    }
                    kernelSum = kernelSum > 255 ? 255 : kernelSum < 0 ? 0 : kernelSum;
                    fImage[(y - GAP) * imageW + (x - GAP)] = kernelSum;

                }
            });

           byte[] res = new byte[imageW * imageH];

            
           Parallel.For(0, imageH, y =>
           {
               for (int x = 0; x < imageW; x++)
               {
                   res[y * imageW + x] = (byte)fImage[y * imageW + x];
               }
           });

            return res;
        }
        public static double [] /***/HC_FILTER_Convolution(double[] fKernel, double[] fRawImage, int imageW, int imageH)
        {
            int KSIZE = (int)Math.Sqrt(fKernel.Length);
            int GAP = KSIZE / 2;

            double [] rawExpanded = (double[])ARRAY_Padding_ALL(fRawImage, imageW, imageH, GAP);

            double [] fRawRes = new double[imageW * imageH];

            int imageNewW = imageW + GAP * 2;
            int imageNewH = imageH + GAP * 2;

            Parallel.For(GAP, imageNewH - GAP, y =>
            {
                for (int x = GAP; x < imageNewW - GAP; x++)
                {
                    double kernelSum = 0;
                    for (int j = -GAP; j <= GAP; j++)
                    {
                        for (int k = -GAP; k <= GAP; k++)
                        {
                            kernelSum += (fKernel[(j + GAP) * KSIZE + k + GAP] * rawExpanded[(y - j) * imageNewW + (x - k)]);
                        }
                    }
                    kernelSum = kernelSum > 255 ? 255 : kernelSum < 0 ? 0 : kernelSum;
                    fRawRes[(y - GAP) * imageW + (x - GAP)] = kernelSum;
                }
            });
            return fRawRes;
        }
        public static byte[] /*****/HC_FILTER_ConvolutionWindow(double[] fKernel, byte[] rawImage, int imageW, int imageH, Rectangle rc)
        {
            double[] fImage = new double[rc.Width*rc.Height];

            int KSIZE = (int)Math.Sqrt(fKernel.Length);
            int GAP = KSIZE / 2;

            //for( int y = rc.Y; y < rc.Y + rc.Height; y++)
            Parallel.For( rc.Y, rc.Y+rc.Height, y=>
            {
                for (int x = rc.X; x < rc.X+rc.Width; x++)
                {
                    double kernelSum = 0;
                    for (int j = -GAP; j <= GAP; j++)
                    {
                        for (int k = -GAP; k <= GAP; k++)
                        {
                            kernelSum += (fKernel[(j + GAP) * KSIZE + k + GAP] * rawImage[(y - j) * imageW + (x - k)]);
                        }
                    }
                    kernelSum = kernelSum > 255 ? 255 : kernelSum < 0 ? 0 : kernelSum;
                    fImage[(y - GAP) * imageW + (x - GAP)] = kernelSum;
                }
            });

            Parallel.For( rc.Y, rc.Y+rc.Height, y =>
            {
                for (int x = rc.X; x < rc.X+rc.Width; x++)
                {
                    rawImage[y * imageW + x] = (byte)fImage[y * imageW + x];
                }
            });

            return rawImage;
        }
        #endregion

        //*****************************************************************************************
        // Gaussian 
        //*****************************************************************************************

        #region GAUSSIAN
        public static double[] HC_FILTER_GenerateGaussianFilter(double fSigma, int nKSize)
        {
            double[] fKernel = new double[nKSize * nKSize];

            int GAP = nKSize / 2;

            //for (int y = -GAP; y <= GAP; y++)
            //{
            //    for (int x = -GAP; x <= GAP; x++)
            //    {
            //        fKernel[(y + GAP) * nKSize + x + GAP] = x;
            //    }
            //}

            double s = 2.0 * fSigma * fSigma;

            double fSum = 0;

            for (int x = -GAP; x <= GAP; x++)
            {
                for (int y = -GAP; y <= GAP; y++)
                {
                    //double r = Math.Sqrt(x * x + y * y);
                    //
                    //fKernel[(y + GAP) * nKSize + x + GAP] = Math.Exp((-((r * r) / s))) / (s * Math.PI);
                    //fSum += fKernel[(y + GAP) * nKSize + x + GAP];

                    fKernel[(y + GAP) * nKSize + x + GAP] = Math.Exp((-((x * x + y * y) / 2*fSigma*fSigma))) / (2 * Math.PI* fSigma * fSigma);
                    fSum += fKernel[(y + GAP) * nKSize + x + GAP];
                }
            }

            for (int y = 0; y < nKSize; y++)
            {
                for (int x = 0; x < nKSize; x++)
                {
                    fKernel[y * nKSize + x] /= fSum;
                }
            }

            return fKernel;
        }
        public static double[] HC_FILTER_GenerateGaussianFilter1(double fSigma, int nKSize)
        {
            double[] fKernel = new double[nKSize * nKSize];

            int GAP = nKSize / 2;

            double fSum = 0;

            for (int x = -GAP; x <= GAP; x++)
            {
                for (int y = -GAP; y <= GAP; y++)
                {

                    fKernel[(y + GAP) * nKSize + x + GAP] = (Math.Exp((-((x * x + y * y) / 2 * fSigma * fSigma))) * (1.0 - ((x * x + y * y) / 2.0 * fSigma * fSigma)) * -(Math.PI * fSigma * fSigma * fSigma * fSigma));

                    //fKernel[(y + GAP) * nKSize + x + GAP] = Math.Exp((-((x * x + y * y) / 2 * fSigma * fSigma))) / (2 * Math.PI * fSigma * fSigma);
                    fSum += fKernel[(y + GAP) * nKSize + x + GAP];
                }
            }

            for (int y = 0; y < nKSize; y++)
            {
                for (int x = 0; x < nKSize; x++)
                {
                    fKernel[y * nKSize + x] /= fSum;
                }
            }

            return fKernel;
        }
        public static double[] HC_FILTER_GenerateLogFilter(int nKSize)
        {
            double[] fKernel = null;
            if (nKSize == 3)
            {
                fKernel = new double []{ -1, -1, -1, -1, 8, -1, -1, -1, -1 };
            }
            else if (nKSize == 5)
            {
                fKernel = new double[]{ 0,0,  1,0,0, 
                                        0,1,  2,1,0, 
                                        1,2,-16,2,1, 
                                        0,1,  2,1,0, 
                                        0,0,  1,0,0};
            }
            else if (nKSize == 7)
            {
                fKernel = new double[] { 0, 0,  1,   1,  1, 0, 0, 
                                         0, 1,  3,   3,  3, 1, 0, 
                                         1, 3,  0, - 7,  0, 3, 1, 
                                         1, 3, -7, -24, -7, 3, 1, 
                                         1, 3,  0, - 7,  0, 3, 1, 
                                         0, 1,  3,   3,  3, 1, 0, 
                                         0, 0,  1,   1,  1, 0, 0 };
            }
            else if (nKSize == 9)
            {
                fKernel = new double[] { 0,0,3,  2,  2,  2,3,0,0,
                                         0,2,3,  5,  5,  5,3,2,0,
                                         3,3,5,  3,  0,  3,5,3,3,
                                         2,5,3,-12,-23,-12,3,5,2,
                                         2,5,0,-23,-40,-23,0,5,2,
                                         2,5,3,-12,-23,-12,3,5,2,
                                         3,3,5,  3,  0,  3,5,3,3,
                                         0,2,3,  5,  5,  5,3,2,0,
                                         0,0,3,  2,  2,  2,3,0,0};
            }
            return fKernel;
        }
        #endregion

        //*****************************************************************************************
        // Prewitt
        //*****************************************************************************************

        #region PREWITT BINARY

        public static double[] HC_Filter_GeneratePrewitt_X(int nKSize)
        {
            double[] fKernel = null;

            if (nKSize == 3)
            {
                fKernel = new double []{ -1, 0, 1, -1, 0, 1, -1, 0, 1 };
            }
            else if (nKSize == 5)
            {
                fKernel = new double[] { 9, 9, 9, 9, 9, 9, 5, 5, 5, 9, -7, -3, 0, -3, -7, -7, -3, -3, -3, -7, -7, -7, -7, -7, -7 };
            }

            return fKernel;
        }
        public static byte[] HC_FILTER_Prewitt_Hor_Bin(byte[] rawImage, int imageW, int imageH)
        {
            double[] kernel = { 1, 1, 1, 0, 0, 0, -1, -1, -1 };

            byte [] res = HC_FILTER_Convolution(kernel, rawImage, imageW, imageH);

            return res;
        }
        public static byte[] HC_FILTER_Prewitt_Ver_Bin(byte[] rawImage, int imageW, int imageH)
        {
            double[] kernel = { -1,0,1,-1,0,1,-1,0,1 };

            byte [] res = HC_FILTER_Convolution(kernel, rawImage, imageW, imageH);

            return res;
        }
        public static byte[] HC_FILTER_Prewitt_Bin(byte[] rawImage, int imageW, int imageH)
        {
            byte[] rawImageHor = HC_FILTER_Prewitt_Hor_Bin(rawImage, imageW, imageH);
            byte[] rawImageVer = HC_FILTER_Prewitt_Ver_Bin(rawImage, imageW, imageH);

            return HC_ARITH_OR(rawImageHor, rawImageVer, imageW, imageH);
        }
        #endregion

        #region PREWITT GRAY
        public static byte[] HC_FILTER_Prewitt_Hor_Gray(byte[] rawImage, int imageW, int imageH)
        {
            double[] kernel = { 1, 1, 1, 0, 0, 0, -1, -1, -1 };

            double[] fImage = HC_CONV_Byte2Double(rawImage);
            fImage = HC_FILTER_Convolution(kernel, fImage, imageW, imageH);

            return HC_CONV_Double2Byte(fImage);
        }
        public static byte[] HC_FILTER_Prewitt_Ver_Gray(byte[] rawImage, int imageW, int imageH)
        {
            double[] kernel = { -1, 0, 1, -1, 0, 1, -1, 0, 1 };

            double[] fImage = HC_CONV_Byte2Double(rawImage);
            fImage = HC_FILTER_Convolution(kernel, fImage, imageW, imageH);

            return HC_CONV_Double2Byte(fImage);
        }
        public static byte[] HC_FILTER_Prewitt_Gray(byte[] rawImage, int imageW, int imageH)
        {
            byte[] rawImageHor = HC_FILTER_Prewitt_Hor_Gray(rawImage, imageW, imageH);
            byte[] rawImageVer = HC_FILTER_Prewitt_Ver_Gray(rawImage, imageW, imageH);

            return HC_ARITH_MID(rawImageHor, rawImageVer, imageW, imageH);
        }
        #endregion
        //*****************************************************************************************
        // Sobel 
        //*****************************************************************************************

        #region SOBEL BINARY

        public static byte[] HC_FILTER_Sobel_Hor_Bin(byte[] rawImage, int imageW, int imageH)
        {
            double[] kernel = { 1, 2, 1, 0, 0, 0, -1, -2, -1 };

            byte[] res = HC_FILTER_Convolution(kernel, rawImage, imageW, imageH);

            return res;
        }
        public static byte[] HC_FILTER_Sobel_Ver_Bin(byte[] rawImage, int imageW, int imageH)
        {
            double[] kernel = { -1, 0, 1, -2, 0,2, -1, 0, 1 };

            byte[] res = HC_FILTER_Convolution(kernel, rawImage, imageW, imageH);

            return res;
        }
        public static byte[] HC_FILTER_Sobel_Bin(byte[] rawImage, int imageW, int imageH)
        {
            byte[] rawImageHor = HC_FILTER_Sobel_Hor_Bin(rawImage, imageW, imageH);
            byte[] rawImageVer = HC_FILTER_Sobel_Ver_Bin(rawImage, imageW, imageH);

            return HC_ARITH_OR(rawImageHor, rawImageVer, imageW, imageH);
        }
        #endregion

        #region SOBEL GRAY
        public static byte[] HC_FILTER_Sobel_Hor_Gray(byte[] rawImage, int imageW, int imageH)
        {
            double[] kernel = { 1, 2, 1, 0, 0, 0, -1, -2, -1 };

            double[] fImage = HC_CONV_Byte2Double(rawImage);
            fImage = HC_FILTER_Convolution(kernel, fImage, imageW, imageH);

            return HC_CONV_Double2Byte(fImage);
        }
        public static byte[] HC_FILTER_Sobel_Ver_Gray(byte[] rawImage, int imageW, int imageH)
        {
            double[] kernel = { -1, 0, 1, -2, 0, 2, -1, 0, 1 };

            double[] fImage = HC_CONV_Byte2Double(rawImage);
            fImage = HC_FILTER_Convolution(kernel, fImage, imageW, imageH);

            return HC_CONV_Double2Byte(fImage);
        }
        public static byte[] HC_FILTER_Sobel_Gray(byte[] rawImage, int imageW, int imageH)
        {
            byte[] rawImageHor = HC_FILTER_Sobel_Hor_Gray(rawImage, imageW, imageH);
            byte[] rawImageVer = HC_FILTER_Sobel_Ver_Gray(rawImage, imageW, imageH);

            return HC_ARITH_MID(rawImageHor, rawImageVer, imageW, imageH);
        }
        #endregion
        
        //*****************************************************************************************
        // Kirsch
        //*****************************************************************************************

        #region Kirsch BINARY
        public static byte[] HC_FILTER_Kirsch_Hor_Bin(byte[] rawImage, int imageW, int imageH)
        {
            double[] kernel = { 5, 5, 5, -3, 0, -3, -3, -3, -3 };

            byte[] res = HC_FILTER_Convolution(kernel, rawImage, imageW, imageH);

            return res;
        }
        public static byte[] HC_FILTER_Kirsch_Ver_Bin(byte[] rawImage, int imageW, int imageH)
        {
            double[] kernel = { -3, -3, 5, -3, 0, 5, -3, -3, 5 };

            byte[] res = HC_FILTER_Convolution(kernel, rawImage, imageW, imageH);

            return res;
        }
        public static byte[] HC_FILTER_Kirsch_Hor_Gray(byte[] rawImage, int imageW, int imageH)
        {
            double[] kernel = { 5, 5, 5, -3, 0, -3, -3, -3, -3 };

            double [] fImage = HC_CONV_Byte2Double(rawImage);
            fImage = HC_FILTER_Convolution(kernel, fImage, imageW, imageH);

            return HC_CONV_Double2Byte(fImage);
        }
        #endregion

        #region Kirsch GRAY
        public static byte[] HC_FILTER_Kirsch_Ver_Gray(byte[] rawImage, int imageW, int imageH)
        {
            double[] kernel = { -3, -3, 5, -3, 0, 5, -3, -3, 5 };

            double[] fImage = HC_CONV_Byte2Double(rawImage);
            fImage = HC_FILTER_Convolution(kernel, fImage, imageW, imageH);

            return HC_CONV_Double2Byte(fImage);
        }
        public static byte[] HC_FILTER_Kirsch_BIN(byte[] rawImage, int imageW, int imageH)
        {
            byte[] rawImageHor = HC_FILTER_Kirsch_Hor_Bin(rawImage, imageW, imageH);
            byte[] rawImageVer = HC_FILTER_Kirsch_Ver_Bin(rawImage, imageW, imageH);

            return HC_ARITH_OR(rawImageHor, rawImageVer, imageW, imageH);
        }
        public static byte[] HC_FILTER_Kirsch_Gray(byte[] rawImage, int imageW, int imageH)
        {
            byte[] rawImageHor = HC_FILTER_Kirsch_Hor_Gray(rawImage, imageW, imageH);
            byte[] rawImageVer = HC_FILTER_Kirsch_Ver_Gray(rawImage, imageW, imageH);

            return HC_ARITH_MID(rawImageHor, rawImageVer, imageW, imageH);
        }
        #endregion

        //*****************************************************************************************
        // Laplacian of Gaussian
        //*****************************************************************************************

        
        public static byte[] HC_FILTER_LofGauss_BIN(byte[] rawImage, int imageW, int imageH, int nKernelSize = 5)
        {
            double[] kernel = HC_FILTER_GenerateLogFilter(nKernelSize);
            
            byte[] res = HC_FILTER_Convolution(kernel, rawImage, imageW, imageH);

            return res;
        }
      

        #region LAPLACIAN of GAUSSIAN GRAY
        public static byte[] HC_FILTER_LofGauss_Gray(byte[] rawImage, int imageW, int imageH, int nKernelSize = 5)
        {
            double[] kernel = HC_FILTER_GenerateLogFilter(nKernelSize);

            double[] fImage = HC_CONV_Byte2Double(rawImage);
            fImage = HC_FILTER_Convolution(kernel, fImage, imageW, imageH);

            return HC_CONV_Double2Byte(fImage);  
        }
        #endregion

        //*****************************************************************************************
        // Sharpen Edge
        //*****************************************************************************************

        #region SHARPEN EDGE BINARY
        private static double[] HC_FILTER_GetFilterSharpen(int nKernelSize = 3)
        {
            double [] k = new double[nKernelSize*nKernelSize];

            int nCenter = nKernelSize /2;

            for (int i = 0; i < k.Length; i++)
            {
                k[i] = -1;
            }
            k[nCenter * nKernelSize + nCenter] = k.Length;

            return k;
        }
        public static byte[] HC_FILTER_Sharpen_BIN(byte [] rawImage, int imageW, int imageH, int nKernelSize = 3)
        {
            double[] kernel = HC_FILTER_GetFilterSharpen(nKernelSize);


            byte[] res = HC_FILTER_Convolution(kernel, rawImage, imageW, imageH);

            return res;
        }
        #endregion

        #region SHARPEN EDGE GRAY
        public static byte[] HC_FILTER_Sharpen_Gray(byte[] rawImage, int imageW, int imageH, int nKernelSize = 3)
        {
            double[] kernel = HC_FILTER_GetFilterSharpen(nKernelSize);

            HC_ARRAY_Dump("c:\\a.txt", kernel, nKernelSize, nKernelSize);

            double[] fImage = HC_CONV_Byte2Double(rawImage);
            fImage  = HC_FILTER_Convolution(kernel, fImage, imageW, imageH);

            return HC_CONV_Double2Byte(fImage); ;
        }
        #endregion

        public static byte[] HC_FILTER_Median(byte[] rawImage, int imageW, int imageH, int nKernelSize)
        {
            int nMedian = nKernelSize*nKernelSize / 2;
            int GAP = nKernelSize / 2;

            byte[] rawRes = new byte[imageW * imageH];

            byte[] rawExpanded = (byte[])ARRAY_Padding_ALL(rawImage, imageW, imageH, GAP);
            int imageNewW = imageW + GAP * 2;
            int imageNewH = imageH + GAP * 2;

            
            Parallel.For(GAP, imageNewH - GAP, y =>
            {
                for (int x = GAP ; x < imageNewW - GAP; x++)
                {
                    double[] arrWindow = new double[nKernelSize*nKernelSize];

                    int index = 0;
                    for (int h = -GAP ; h <= GAP; h++)
                    {
                        for (int w = -GAP; w <= GAP; w++, index++)
                        {
                            arrWindow[index] = rawExpanded[(y+h) * imageNewW + (x+w)];
                        }
                    }
                    Array.Sort(arrWindow);
                    rawRes[(y - GAP) * imageW + (x - GAP)] = Convert.ToByte(arrWindow[nMedian]);
                }
            });
            return rawRes;
        }
        public static byte[] HC_FILTER_MaskMedian(byte[] mask, byte[] rawImage, int imageW, int imageH, int nKernelSize, int nonTargetValue)
        {
            int nMedian = nKernelSize * nKernelSize / 2;
            int GAP = nKernelSize / 2;

            int nCenterX = imageW / 2;
            int nCenterY = imageH / 2;

            //int nRadius = Convert.ToInt32(Math.Sqrt((imageW * imageW) + (imageH * imageH)) / 2.0);

            byte[] newRaw = new byte[imageW * imageH];

            for (int y = 0; y < imageH; y++)
            {
                for (int x = 0; x < imageW; x++)
                {
                    if (Math.Sqrt((nCenterX - x)*(nCenterX - x)) + ((nCenterY - y)*(nCenterY - y)) > 140)
                    {
                        newRaw[y * imageW + x] = (byte)nonTargetValue;
                    }
                    else
                    {
                        if (mask[y * imageW + x] == 0)
                        {
                            double[] arrWindow = new double[nKernelSize * nKernelSize];
                            int index = 0;
                            for (int wy = y - GAP; wy < y + GAP; wy++)
                            {
                                for (int wx = x - GAP; wx < x + GAP; wx++)
                                {

                                    arrWindow[index++] = rawImage[wy * imageW + wx];


                                }
                            }
                            Array.Sort(arrWindow);
                            newRaw[y * imageW + x] = Convert.ToByte(arrWindow[nMedian]);
                        }
                        else
                        {
                            newRaw[y * imageW + x] = rawImage[y * imageW + x];
                        }

                    }
                }
            }


            return newRaw;
        }

//        public static void CreateSobelKernel(int n, ref float[][] Kx, ref float[][] Ky)
//        {
//            int side = n * 2 + 3;
//            int halfSide = side / 2;
//            for (int i = 0; i < side; i++)
//            {
//                int k = (i <= halfSide) ? (halfSide + i) : (side + halfSide - i - 1);
//                for (int j = 0; j < side; j++)
//                {
//                    if (j < halfSide)
//                        Kx[i][j] = Ky[j][i] = j - k;
//                    else if (j > halfSide)
//                        Kx[i][j] = Ky[j][i] = k - (side - j - 1);
//                    else
//                        Kx[i][j] = Ky[j][i] = 0;
//                }
//            }
//        }

        //*****************************************************************************************
        // Mopology
        //*****************************************************************************************

        #region MOPOLOGY  
        public static byte[] HC_FILTER_Dilation(byte[] rawImage, int imageW, int imageH, int nKernelSize)
        {
            int nGap = nKernelSize / 2;

            byte[] newRes = new byte[imageW * imageH];

            byte[] rawExpanded = (byte[])ARRAY_Padding_ALL(rawImage, imageW, imageH, nGap);
            int imageNewW = imageW + nGap * 2;
            int imageNewH = imageH + nGap * 2;

            
            //for (int y = nGap; y < imageNewH - nGap; y++)
            Parallel.For(nGap, imageNewH - nGap, y =>
            {
                for (int x = nGap; x < imageNewW - nGap; x++)
                {
                    byte max = 0;
                    for (int h = -nGap; h <= nGap; h++)
                    {
                        for (int w = -nGap; w <= nGap; w++)
                        {
                            if (rawExpanded[(y + h) * imageNewW + (x + w)] > max)
                            {
                                max = rawExpanded[(y + h) * imageNewW + (x + w)];
                            }
                        }
                    }
                    newRes[(y - nGap) * imageW + (x - nGap)] = max;
                }
            });

            return newRes;
        }
        public static byte[] HC_FILTER_Erosion(byte[] rawImage, int imageW, int imageH, int nKernelSize)
        {
            int nGap = nKernelSize / 2;

            byte[] newRes = new byte[imageW * imageH];

            byte[] rawExpanded = (byte[])ARRAY_Padding_ALL(rawImage, imageW, imageH, nGap);
            int imageNewW = imageW + nGap * 2;
            int imageNewH = imageH + nGap * 2;

            for (int y = nGap; y < imageNewH - nGap; y++)
            {
                for (int x = nGap; x < imageNewW - nGap; x++)
                {
                    byte min = byte.MaxValue;
                    for (int h = - nGap; h <=  nGap; h++)
                    {
                        for (int w = - nGap; w <=  nGap; w++)
                        {
                            if (rawExpanded[(y+h) * imageNewW + (x+w)] < min)
                            {
                                min = rawExpanded[(y+h) * imageNewW+ (x+w)];
                            }
                        }
                    }
                    newRes[(y-nGap) * imageW + (x-nGap)] = min;
                }
            }

            return newRes;
        }
        public static byte[] HC_FILTER_Open(byte[] rawImage, int imageW, int imageH, int nKernelSize)
        {
            byte[] newRes = null;
            Parallel.Invoke(() =>{newRes = HC_FILTER_Erosion(rawImage, imageW, imageH, nKernelSize);});
            Parallel.Invoke(() =>{newRes = HC_FILTER_Dilation(newRes, imageW, imageH, nKernelSize);});
            return newRes;
        }
        public static byte[] HC_FILTER_Close(byte[] rawImage, int imageW, int imageH, int nKernelSize)
        {
            byte[] newRes = null;
            Parallel.Invoke(() =>{newRes = HC_FILTER_Dilation(rawImage, imageW, imageH, nKernelSize);});
            Parallel.Invoke(() =>{newRes = HC_FILTER_Erosion(newRes, imageW, imageH, nKernelSize); });
            return newRes;
        }
        #endregion 
    }
}


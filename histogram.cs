using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Drawing;
using System.Drawing.Imaging;
using DispObject;

namespace ReinlessLib
{
    public static partial class CReinlessLib
    {

        #region HISTOGRAM 

        public static int[] /**/HC_HISTO_GetHistogram(byte[] rawImage, int imageW, int imageH, int nNormalizationValue)
        {
            int[] nHistogram = new int[256];

            for (int y = 0; y < imageH; y++)
            {
                for (int x = 0; x < imageW; x++)
                {
                    int nValue = rawImage[y * imageW + x];
                    nHistogram[nValue]++;
                }
            }
            if( nNormalizationValue > 0 )
                HC_HISTO_Normalization(nHistogram, nHistogram.Length, nNormalizationValue);

            return nHistogram;
        }
        public static int[] HC_HISTO_GetHistogram(byte[] array, int nNormalization)
        {
            int[] nHistogram = new int[256];
            for (int i = 0; i < array.Length; i++)
            {
                nHistogram[array[i]]++;
            }
            if (nNormalization > 0)
            {
                HC_HISTO_Normalization(nHistogram, array.Length, nNormalization);
            }
            return nHistogram;
        }

        public static void /***/HC_HISTO_Normalization(int[] nArray, int nLength, int nNormalizationValue)
        {
            try
            {
                // Fucking Normalization
                int Max = 1;

                for (int nIndex = 0; nIndex < nArray.Length; nIndex++)
                {
                    double fValue = nArray[nIndex];
                    if (fValue > Max) Max = (int)fValue;
                }

                for (int nIndex = 0; nIndex < nArray.Length; nIndex++)
                {
                    double fValue = nArray[nIndex] * nNormalizationValue;
                    nArray[nIndex] = (int)fValue;

                    // i want to avoid divide by zero
                    if (nArray[nIndex] != 0) nArray[nIndex] /= Max;
                }
            }
            catch (Exception ex)
            {
                Console.Write(ex.ToString());
            }

        }
        public static void /***/HC_HISTO_Normalization(double[] fArray, int nLength, int nNormalizationValue)
        {
            try
            {
                // Fucking Normalization
                double Max = 1;

                for (int nIndex = 0; nIndex < fArray.Length; nIndex++)
                {
                    double fValue = fArray[nIndex];
                    if (fValue > Max) Max = fValue;
                }

                for (int nIndex = 0; nIndex < fArray.Length; nIndex++)
                {
                    double fValue = fArray[nIndex] * nNormalizationValue;
                    fArray[nIndex] = fValue;

                    // i want to avoid divide by zero
                    if (fArray[nIndex] != 0) fArray[nIndex] /= Max;
                }
            }
            catch (Exception ex)
            {
                Console.Write(ex.ToString());
            }

        }

        #endregion

        #region Image Statistics

        public static int /******/HC_HISTO_GetMax(byte[] rawImage, int imageW, int imageH)
        {
            int[] Histogram = HC_HISTO_GetHistogram(rawImage, imageW, imageH, 255);
            return Histogram.Max();
        }
        public static int /******/HC_HISTO_GetMin(byte[] rawImage, int imageW, int imageH)
        {
            int[] Histogram = HC_HISTO_GetHistogram(rawImage, imageW, imageH, 255);
            return Histogram.Min();
        }
        public static double /***/HC_HISTO_GetAvg(byte[] rawImage, int imageW, int imageH)
        {
            int[] Histogram = HC_HISTO_GetHistogram(rawImage, imageW, imageH, 255);
            return Histogram.Average();
        }
        public static double /***/HC_HISTO_GetVar(byte[] rawImage, int imageW, int imageH)
        {
            int[] Histogram = HC_HISTO_GetHistogram(rawImage, imageW, imageH, 255);

            double fMin = Histogram.Min();
            double fMax = Histogram.Max();

            return fMax - fMin;
        }
        public static double /***/HC_HISTO_GetStatics(byte[] rawImage, int imageW, int imageH, ref double min, ref double max, ref double var)
        {
            int HCount = imageH / 10;

            int[] Line = new int[imageW * HCount];

            decimal f = 0;
            for (int y = 0; y < imageH; y += HCount)
            {
                for (int i = 0; i < imageW; i++)
                {
                    Line[i] = rawImage[y * imageW + i];
                    f += Line[i];
                }
            }

            f /= Line.Length;
            min = Line.Min();
            max = Line.Max();
            var = max - min;
            return (double)f;
        }

        #endregion

        #region PROJECTION related

        public static int[] /****/HC_Histo_GetProjectionAll(byte[] rawImage, int imageW, int imageH)
        {
            int[] projH = HC_HISTO_GetProjectionH(rawImage, imageW, imageH);
            int[] projV = HC_HISTO_GetProjectionV(rawImage, imageW, imageH);

            int[] proj = new int[imageW + imageH];

            Array.Copy(projH, 0, proj, 0, imageW);
            Array.Copy(projV, 0, proj, imageW, imageH);

            return proj;
        }
        public static int[] /****/HC_HISTO_GetProjectionH(byte [] rawImage, int imageW, int imageH)
         {
            int [] arrProjectionH = new int[imageW];

              int nSum = 0;
              for( int x = 0 ; x < imageW; x++)
              {
                  nSum = 0;
                 for( int y = 0; y < imageH; y++) {nSum += rawImage[y*imageW+x];}
                 arrProjectionH[x] = nSum;
              }
              return arrProjectionH;
         }
        public static int[] /****/HC_HISTO_GetProjectionV(byte[] rawImage, int imageW, int imageH)
        {
            int[] arrProjectionV = new int[imageH];

              int nSum = 0;
              for( int y = 0; y < imageH; y++)
              {
                  nSum = 0;
                  for( int x = 0; x < imageW; x++){nSum += rawImage[y*imageW+x];}
                  arrProjectionV[y] = nSum;
              }
            return arrProjectionV;
         }

        #endregion

        // Get Histogram Distribution Pos [head, center, tail]
        public static void /*****/HC_HISTO_GetDataRange_ALL(int[] nArr, ref int head, ref int mid, ref int tail)
        {
            for (int i = 0; i < nArr.Length; i++){if (nArr[i] != 0){head = i;break;}}
            for (int i = nArr.Length - 1; i >= 0; i--){if (nArr[i] != 0) { tail = i; break; }}
            
            mid = (int)(head + (tail - head)/2.0);
        }
        public static void /*****/HC_HISTO_GetDataRange_HeadTail(int[] nArr, ref int head, ref int tail)
        {
            for (int i = 0; i < nArr.Length; i++) { if (nArr[i] != 0) { head = i; break; } }
            for (int i = nArr.Length - 1; i >= 0; i--) { if (nArr[i] != 0) { tail = i; break; } }
        }

        public static double /***/HC_HISTO_GetIntersection(int[] QueryGrayHisto, int[] TargetGrayHisto)
    {
        int x = 0;
        double Buffer = 0;
        float SumOfModel = 0;

        //*********************************************************************************************
        //		The Histogram of the intersection
        //*********************************************************************************************
        for (x = 0; x < 256; x++)
        {
            if (QueryGrayHisto[x] < TargetGrayHisto[x])
            {
                Buffer += QueryGrayHisto[x];
            }
            else if (QueryGrayHisto[x] > TargetGrayHisto[x])
            {
                Buffer += TargetGrayHisto[x];
            }
            else if (QueryGrayHisto[x] == TargetGrayHisto[x])
            {
                Buffer += QueryGrayHisto[x];
            }
            SumOfModel += TargetGrayHisto[x];
        }
        //*********************************************************************************************
        // To obtain a fractional match value between 0 to 1
        // the intersection is normalized by the number of pixels in the model histogram;
        //*********************************************************************************************

        Buffer = Buffer / SumOfModel;

        return Buffer;
    }
        public static double /***/HC_FOCUS_GetFocusValue(byte[] rawImage, int imageW, int imageH)
        {
            double sum = 0;
            for (int y = 0; y < imageH; y++)
            {
                for (int x = 0; x < imageW - 1; x++)
                {
                    sum += Math.Abs(rawImage[y * imageW + x] - rawImage[y * imageW + x + 1]);
                }
            }

            sum /= imageW * imageH;
            return sum;
        }

        
        
        public static double[] HC_HISTO_GetGradientVectorImage(byte[] rawImage, int imageW, int imageH)
        {
            double[] fImage = new double[imageW * imageH];

            Parallel.For(1, imageH - 1, y =>
            {
                for (int x = 1; x < imageW - 1; x++)
                {
                    double dx = (rawImage[y * imageW + x + 1] - rawImage[y * imageW + x - 1]);
                    double dy = (rawImage[(y + 1) * imageW + x] - rawImage[(y - 1) * imageW + x]);

                    double fValue = Math.Atan(dy / dx) * 180.0 / Math.PI;

                    if (double.IsNaN(fValue))
                        fValue = 0;

                    if (fValue < 0) fValue = 360 + fValue;
                    fImage[y * imageW + x] = fValue;
                }
            });
            return fImage;
        }

        // Get Gradient Angular Histogram
        public static int[] HC_HISTO_GetGradientVector(byte[] rawImage, int imageW, int imageH)
        {
            int [] GradientVector = new int[360];

            //for (int y = 1; y < imageH - 1; y++)
            Parallel.For(1, imageH - 1, y =>
            {
                for (int x = 1; x < imageW - 1; x++)
                {
                    double dx = (rawImage[y * imageW + x + 1] - rawImage[y * imageW + x - 1]);
                    double dy = (rawImage[(y + 1) * imageW + x] - rawImage[(y - 1) * imageW + x]);

                    double fValue = Math.Atan(dy / dx) * 180.0 / Math.PI;

                    if (double.IsNaN(fValue))
                        fValue = 0;

                    if (fValue < 0) fValue = 360 + fValue;

                    GradientVector[(int)fValue]++;
                }
            });
            return GradientVector;

        }

        public static List<DPoint> HC_HISTO_GetGradientPixels(byte[] rawImage, int imageW, int imageH)
        {
            List<DPoint> ptList = new List<DPoint>();

            Parallel.For(1, imageH - 1, y =>
            {
                for (int x = 1; x < imageW - 1; x++)
                {
                    double dx = (rawImage[y * imageW + x + 1] - rawImage[y * imageW + x - 1]);
                    double dy = (rawImage[(y + 1) * imageW + x] - rawImage[(y - 1) * imageW + x]);

                    double fValue = Math.Atan(dy / dx) * 180.0 / Math.PI;

                    if (double.IsNaN(fValue))
                        fValue = 0;

                    if (fValue > 40 && fValue < 60)
                    {
                        ptList.Add(new DPoint(x, y));
                    }
                     
                }
            });

            return ptList;
        }

       
        public static int [] HC_Histo_GetAutoCorrelogram(byte[] rawImage, int imageW, int imageH)
        {
            int[] histo1D = new int[256];
            int[] histo2D = new int[256];

            for (int y = 2; y < imageH - 2; y++)
            {
                for (int x = 2; x < imageW - 2; x++)
                {
                    int nCurrent = rawImage[y * imageW + x];

                    int nMatchD2 = 0;
                    int nMatchD1 = 0;

                    for (int xx = x - 2; xx < x + 2; xx++)
                    {
                        if (rawImage[(y - 2) * imageW + xx] == nCurrent)nMatchD2++;
                        if (rawImage[(y + 2) * imageW + xx] == nCurrent) nMatchD2++;
                    }
                    for (int yy = y - 1; yy < y + 1; yy++)
                    {
                        if (rawImage[yy * imageW + x - 2] == nCurrent) nMatchD2++;
                        if (rawImage[yy * imageW + x + 2] == nCurrent) nMatchD2++;
                    }

                    histo2D[nCurrent] += nMatchD2;

                    for (int yy = y - 1; yy < y + 1; yy++)
                    {
                        for (int xx = x - 1; xx < x + 1; xx++)
                        {
                            if (yy == 0 && xx == 0) continue;

                            if (nCurrent == rawImage[yy * imageW + xx]) nMatchD1++;
                        }
                    }

                    histo1D[nCurrent] += nMatchD1;
                }
            }
            int[] returnHisto = new int[512];

            Array.Copy(histo1D, 0, returnHisto, 0, 256);
            Array.Copy(histo2D, 0, returnHisto, 256, 256);

            return returnHisto;
        }

        
    }
}

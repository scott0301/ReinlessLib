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
        #region  Auto Threshold Approach

        public static byte[]/****/HC_THR_Otsu(byte[] rawImage, int imageW, int imageH)
        {
            int[] nHistogram = HC_HISTO_GetHistogram(rawImage, imageW, imageH, 255);

            int nThreshold = HC_THR_Otsu(nHistogram, imageW, imageH);

            HC_THR_Binarization_Single(rawImage, imageW, imageH, nThreshold);
            return rawImage;
        }
        public static int /******/HC_THR_Otsu(int[] nHistogram, int imageW, int imageH)
        {
            // 스레스홀딩 관련 변수
            double[] betweenVariance = new double[256];     // 클래스간 분산 값들 저장 배열
            double fMaxValue = 0;                            // 최대치
            int nThreshold = 0;                                  // 임계점

            double weightB, meanB;                          // 배경
            double weightF, meanF;                          // 전경
            double temp;

            for (int nBin = 0; nBin < 256; nBin++)
            {
                weightB = 0; meanB = 0;
                weightF = 0; meanF = 0;

                // 배경 : 가충지, 평균값 계산
                for (int nIndex = 0; nIndex <= nBin; nIndex++)
                {
                    weightB = weightB + nHistogram[nIndex];
                    meanB = meanB + (nIndex * nHistogram[nIndex]);
                }

                // Divide By Zero
                temp = weightB;
                if (temp == 0)
                    temp = 0.00001;
                else if (temp == 1)
                    temp = 0.99999;

                weightB = weightB / (imageW * imageH);
                meanB = meanB / temp;

                // 전경 : 가중치, 평균값 계산
                for (int i = nBin + 1; i < 256; i++)
                {
                    weightF = weightF + nHistogram[i];
                    meanF = meanF + (i * nHistogram[i]);
                }

                // Divide By Zero
                temp = weightF;
                if (temp == 0)
                    temp = 0.00001;
                else if (temp == 1)
                    temp = 0.99999;

                weightF = weightF / (imageW * imageH);
                meanF = meanF / temp;

                // 클래스간 분산 계산
                betweenVariance[nBin] = weightB * weightF * (meanB - meanF) * (meanB - meanF);
            }

            // 임계값 찾기
            fMaxValue = betweenVariance[0];
            for (int i = 0; i < 256; i++)
            {
                if (fMaxValue < betweenVariance[i])
                {
                    fMaxValue = betweenVariance[i];
                    nThreshold = i;
                }
            }
            return nThreshold;
        }
        public static int /******/HC_THR_Huang(int[] nHistogram)
        {
            // Implements Huang's fuzzy thresholding method 
            // Uses Shannon's entropy function (one can also use Yager's entropy function) 
            // Huang L.-K. and Wang M.-J.J. (1995) "Image Thresholding by Minimizing  
            // the Measures of Fuzziness" Pattern Recognition, 28(1): 41-51
            // M. Emre Celebi  06.15.2007

            int threshold = -1;
            int ih, it;
            int sum_pix;
            int num_pix;
            double term;
            double ent;  // entropy 
            double min_ent; // min entropy 
            double mu_x;

            /* Determine the first non-zero bin */
            int first_bin = 0;
            for (ih = 0; ih < 256; ih++)
            {
                if (nHistogram[ih] != 0)
                {
                    first_bin = ih;
                    break;
                }
            }

            /* Determine the last non-zero bin */
            int last_bin = 255;
            for (ih = 255; ih >= first_bin; ih--)
            {
                if (nHistogram[ih] != 0)
                {
                    last_bin = ih;
                    break;
                }
            }
            term = 1.0 / (double)(last_bin - first_bin);
            double[] mu_0 = new double[256];
            sum_pix = num_pix = 0;
            for (ih = first_bin; ih < 256; ih++)
            {
                sum_pix += ih * nHistogram[ih];
                num_pix += nHistogram[ih];
                /* NUM_PIX cannot be zero ! */
                mu_0[ih] = sum_pix / (double)num_pix;
            }

            double[] mu_1 = new double[256];
            sum_pix = num_pix = 0;
            for (ih = last_bin; ih > 0; ih--)
            {
                sum_pix += ih * nHistogram[ih];
                num_pix += nHistogram[ih];
                /* NUM_PIX cannot be zero ! */
                mu_1[ih - 1] = sum_pix / (double)num_pix;
            }

            /* Determine the threshold that minimizes the fuzzy entropy */
            threshold = -1;
            min_ent = Double.MaxValue;
            for (it = 0; it < 256; it++)
            {
                ent = 0.0;
                for (ih = 0; ih <= it; ih++)
                {
                    /* Equation (4) in Ref. 1 */
                    mu_x = 1.0 / (1.0 + term * Math.Abs(ih - mu_0[it]));
                    if (!((mu_x < 1e-06) || (mu_x > 0.999999)))
                    {
                        /* Equation (6) & (8) in Ref. 1 */
                        ent += nHistogram[ih] * (-mu_x * Math.Log(mu_x) - (1.0 - mu_x) * Math.Log(1.0 - mu_x));
                    }
                }

                for (ih = it + 1; ih < 256; ih++)
                {
                    /* Equation (4) in Ref. 1 */
                    mu_x = 1.0 / (1.0 + term * Math.Abs(ih - mu_1[it]));
                    if (!((mu_x < 1e-06) || (mu_x > 0.999999)))
                    {
                        /* Equation (6) & (8) in Ref. 1 */
                        ent += nHistogram[ih] * (-mu_x * Math.Log(mu_x) - (1.0 - mu_x) * Math.Log(1.0 - mu_x));
                    }
                }
                /* No need to divide by NUM_ROWS * NUM_COLS * LOG(2) ! */
                if (ent < min_ent)
                {
                    min_ent = ent;
                    threshold = it;
                }
            }
            return threshold;

        }
        public static int /******/HC_THR_Yen(int[] nHistogram)
        {
            // Implements Yen  thresholding method
            // 1) Yen J.C., Chang F.J., and Chang S. (1995) "A New Criterion 
            //    for Automatic Multilevel Thresholding" IEEE Trans. on Image 
            //    Processing, 4(3): 370-378
            // 2) Sezgin M. and Sankur B. (2004) "Survey over Image Thresholding 
            //    Techniques and Quantitative Performance Evaluation" Journal of 
            //    Electronic Imaging, 13(1): 146-165
            //    http://citeseer.ist.psu.edu/sezgin04survey.html
            //
            // M. Emre Celebi
            // 06.15.2007
            // Ported to ImageJ plugin by G.Landini from E Celebi's fourier_0.8 routines
            int threshold;
            int ih, it;
            double crit;
            double max_crit;
            double[] norm_histo = new double[256]; /* normalized histogram */
            double[] P1 = new double[256]; /* cumulative normalized histogram */
            double[] P1_sq = new double[256];
            double[] P2_sq = new double[256];

            int total = 0;
            for (ih = 0; ih < 256; ih++)
                total += nHistogram[ih];

            for (ih = 0; ih < 256; ih++)
                norm_histo[ih] = (double)nHistogram[ih] / total;

            P1[0] = norm_histo[0];
            for (ih = 1; ih < 256; ih++)
                P1[ih] = P1[ih - 1] + norm_histo[ih];

            P1_sq[0] = norm_histo[0] * norm_histo[0];
            for (ih = 1; ih < 256; ih++)
                P1_sq[ih] = P1_sq[ih - 1] + norm_histo[ih] * norm_histo[ih];

            P2_sq[255] = 0.0;
            for (ih = 254; ih >= 0; ih--)
                P2_sq[ih] = P2_sq[ih + 1] + norm_histo[ih + 1] * norm_histo[ih + 1];

            /* Find the threshold that maximizes the criterion */
            threshold = -1;
            max_crit = Double.MinValue;
            for (it = 0; it < 256; it++)
            {
                crit = -1.0 * ((P1_sq[it] * P2_sq[it]) > 0.0 ? Math.Log(P1_sq[it] * P2_sq[it]) : 0.0) + 2 * ((P1[it] * (1.0 - P1[it])) > 0.0 ? Math.Log(P1[it] * (1.0 - P1[it])) : 0.0);
                if (crit > max_crit)
                {
                    max_crit = crit;
                    threshold = it;
                }
            }
            return threshold;
        }        

        #endregion

        #region Binarization 

        public static byte[] /***/HC_THR_Binarization_Single(byte[] rawImage, int imageW, int imageH, int nThreshold)
        {
            byte[] rawRes = new byte[imageW * imageH];
            for (int i = 0; i < imageW * imageH; i++)
            {
                rawRes[i] = rawImage[i] > nThreshold ? (byte)255 : (byte)0;
            }
            return rawRes;
        }
        public static byte[] /***/HC_THR_Binarization_Dual(byte[] rawImage, int imageW, int imageH, int nThrLow, int nThrHigh)
        {
            byte[] rawRes = new byte[imageW * imageH];

            for (int i = 0; i < imageW * imageH; i++)
            {
                rawRes[i] = rawImage[i] > nThrHigh ? (byte)0 : rawImage[i] > nThrLow ? (byte)255 : (byte)0;
            }
            return rawRes;
        }
        public static byte[] /***/HC_THR_MaximumEntropy(byte[] rawImage, int imageW, int imageH, int start, int end)
        {
            int[] nHistogram = HC_HISTO_GetHistogram(rawImage, imageW, imageH, 0);

            double[] fHistoNorm = new double[256];
            double norm = imageW*imageH;

            Parallel.For(0,  256, i => {fHistoNorm[i] = nHistogram[i] / norm; });

            double[] entropy = new double[256];

            double sum_prob_1k = 0;
            double sum_prob_kl = 0;
            double sum_prob_in_1k = 0;
            double sum_prob_ln_kl = 0;

            //*************************************************************************************
            // Calculate Entropy

            for (int k = start; k < end; k++)
            {
                sum_prob_1k = 0;
                sum_prob_kl = 0;
                sum_prob_in_1k = 0;
                sum_prob_ln_kl = 0;

                for (int i = 1; i < k; i++)
                {
                    sum_prob_1k += fHistoNorm[i];
                    if (fHistoNorm[i] != 0) { sum_prob_in_1k += fHistoNorm[i] * Math.Log(fHistoNorm[i]); }
                }
                for (int i = k; i < end; i++)
                {
                    sum_prob_kl += fHistoNorm[i];
                    if (fHistoNorm[i] != 0) { sum_prob_ln_kl += (fHistoNorm[i] *Math.Log(fHistoNorm[i])); }
                }

                entropy[k] = Math.Log(sum_prob_1k) + Math.Log(sum_prob_kl) - (sum_prob_in_1k/sum_prob_1k) - (sum_prob_ln_kl / sum_prob_kl);

                if( entropy[k] < 0 ) { entropy[k] = 0; }
            }

            //****************************************************************************************
            // Generate Entropy Look up table

            double fMaxEntropy = entropy.Max();
            int nMaxPos = Array.IndexOf(entropy, fMaxEntropy);

            byte[] LUT = new byte[256];

            for (int i = 0; i < 256; i++)
            {
                /***/if (i <= nMaxPos) LUT[i] = 0;
                else if (i > nMaxPos) LUT[i] = 255;
            }

            byte [] newRaw = new byte[imageW*imageH];

            for (int i = 0; i < newRaw.Length; i++) {newRaw[i] = LUT[rawImage[i]];}

            return newRaw;
        }
        #endregion

        // extract arbitrary color = brightness gravity : aussume that we have segmented object which has only one color 170105
        public static PointF /***/HC_THR_OBJECT_GRAVITY(byte[] rawImage, int imageW, int imageH, int Target)
        {
            float cx = 0;
            float cy = 0;
            int nCount = 0;

            for (int y = 0; y < imageH; y++)
            {
                for (int x = 0; x < imageW; x++)
                {
                    if (rawImage[y * imageW + x] > Target)
                    {
                        cx += x;
                        cy += y;
                        nCount++;
                    }
                }
            }

            cx = (int)(cx / nCount);
            cy = (int)(cy / nCount);

            return new Point((int)cx, (int)cy);
        }

        // Extract Inner Blob (black) with out labeling 170106 
        public static byte[] HC_THR_GenerateInnerBlob(byte[] rawImage, int imageW, int imageH)
        {
            byte[] rawBlob = new byte[rawImage.Length];

            for (int x = 0; x < imageW; x++)
            {
                bool bOn = false;

                for (int y = 1; y < imageH; y++)
                {
                    byte prev = rawImage[(y - 1) * imageW + x];
                    byte curr = rawImage[y * imageW + x];
                    // check color changing
                    if (prev != 0 && curr == 0 && bOn == false) 
                    {
                        rawBlob[y * imageW + x] = 0xFF;
                        bOn = true;
                    }
                    else if (prev == 0 && curr == 0 && bOn == true)
                    {
                        rawBlob[y * imageW + x] = 0xFF;
                    }
                    else if (prev == 0 && curr != 0)
                    {
                        bOn = false;
                    }
                }
                for (int y = imageH - 1; y >= 0; y--)
                {
                    if (rawBlob[y * imageW + x] != 0)
                    {
                        rawBlob[y * imageW + x] = 0;
                    }
                    else
                    {
                        y = 0;
                    }
                }// rollback loop
            }// loop x
            return rawBlob;
        }

        // get blob region rectangle by projection analysis 170106
        public static RectangleF HC_THR_GetObjectRect(byte[] rawImage, int imageW, int imageH)
        {
            int[] projH = CReinlessLib.HC_HISTO_GetProjectionH(rawImage, imageW, imageH);
            int[] projV = CReinlessLib.HC_HISTO_GetProjectionV(rawImage, imageW, imageH);

            int H_head = 0; int H_tail = 0;
            int V_head = 0; int V_tail = 0;

            CReinlessLib.HC_HISTO_GetDataRange_HeadTail(projH, ref H_head, ref H_tail);
            CReinlessLib.HC_HISTO_GetDataRange_HeadTail(projV, ref V_head, ref V_tail);

            PointF P1 = new PointF(H_head, V_head);
            PointF P2 = new PointF(H_tail, V_tail);

            RectangleF rc = new RectangleF(P1.X, P1.Y, P2.X - P1.X, P2.Y - P1.Y);

            return rc;
        }
        /// <summary>
        /// 추출된 Rectangle List에서 Overlapped Rectangle Ellimination 
        /// </summary>
        /// <param name="list"></param>
        /// <param name="nInflateSize"></param>
        /// <returns></returns>
        public static List<Rectangle> HC_BLOB_GetMergedList(List<Rectangle> list, int nInflateSize)
        {
            float MERGE_INFLATE = 4;
            float MERGE_HALF = Convert.ToSingle(MERGE_INFLATE / 2.0);

            bool bOverlapped = true;

            while (bOverlapped == true)
            {
                for (int i = 0; i < list.Count; i++)
                {
                    for (int j = 1; j < list.Count; j++)
                    {
                        if (list.ElementAt(i) == list.ElementAt(j)) continue;

                        RectangleF rcOrg = new RectangleF(list.ElementAt(i).X, list.ElementAt(i).Y, list.ElementAt(i).Width, list.ElementAt(i).Height);
                        RectangleF rcTar = new RectangleF(list.ElementAt(j).X, list.ElementAt(j).Y, list.ElementAt(j).Width, list.ElementAt(j).Height);

                        rcOrg.Offset(-MERGE_HALF, -MERGE_HALF); rcOrg.Inflate(MERGE_INFLATE, MERGE_INFLATE);
                        rcTar.Offset(-MERGE_HALF, -MERGE_HALF); rcTar.Inflate(MERGE_INFLATE, MERGE_INFLATE);

                        //if (list.ElementAt(i).IntersectsWith(list.ElementAt(j)) == true)
                        if (rcOrg.IntersectsWith(rcTar))
                        {
                            list[i] = Rectangle.Union(list.ElementAt(i), list.ElementAt(j));
                            list.RemoveAt(j);
                            i = list.Count;
                            break;
                        }
                    }
                }

                int nDupplicate = 0;

                for (int i = 0; i < list.Count; i++)
                {
                    for (int j = 1; j < list.Count; j++)
                    {
                        if (list.ElementAt(i) == list.ElementAt(j)) continue;
                        if (list.ElementAt(i).IntersectsWith(list.ElementAt(j)))
                        {
                            nDupplicate++;
                        }
                    }
                }

                if (nDupplicate == 0)
                {
                    bOverlapped = false;
                }
            }

            return list.ToList();
        }
    }
}

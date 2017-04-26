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

        #region IMAGE SAVE
        public static byte[] /***/LoadImage(string strPath, ref int imageW, ref int imageH)
        {
            Bitmap bmp = (Bitmap)Image.FromFile(strPath);

            byte[] rawImage = HC_CONV_Bmp2Byte(bmp, ref imageW, ref imageH);
            return rawImage;

        }
        public static void /*****/SaveImage(byte[] rawImage, int imageW, int imageH, string strPath)
        {
            Bitmap bmp = (Bitmap)HC_CONV_Byte2Bmp(rawImage, imageW, imageH);

            bmp.Save(strPath);
        }
        public static void /*****/SaveImage(double[] fImage, int imageW, int imageH, string strPath)
        {
            byte[] rawImage = HC_CONV_Double2Byte(fImage);

            Bitmap bmp = (Bitmap)HC_CONV_Byte2Bmp(rawImage, imageW, imageH);
            bmp.Save(strPath);
        }
        #endregion

   
      
        public static void HC_ARRAY_Dump(string strPath, double[] fArray, int w, int h)
        {
            string strBody = string.Empty;


            if (fArray.Length < w * h) return;

            for (int y = 0; y < h; y++)
            {
                for (int x = 0; x < w; x++)
                {

                    strBody += fArray[y * w + x].ToString() + ",";
                }
                strBody += System.Environment.NewLine;
            }

            try
            {

                System.IO.File.WriteAllText(strPath, strBody);
            }
            catch { }
        }
        public static void HC_ARRAY_Dump(string strPath, int[] nArray)
        {

            string strBody = string.Empty;

            for (int i = 0; i < nArray.Length; i++)
            {
                strBody += nArray[i].ToString() + ",";
            }

            try
            {
                System.IO.File.WriteAllText(strPath, strBody);
            }
            catch { }
        }
        public static void HC_ARRAY_Dump(string strPath, byte[] byteArray, int w, int h)
        {
            string strBody = string.Empty;

            for (int y = 0; y < h; y++)
            {
                for (int x = 0; x < w; x++)
                {
                    strBody += byteArray[y * w + x].ToString() + ",";
                }
                strBody += System.Environment.NewLine;
            }
            System.IO.File.WriteAllText(strPath, strBody);
        }


        public static int GetMaxElementPosition(double[] array)
        {
            double fMax = array.Max();
            int nIndex = Array.IndexOf(array, fMax);

            return nIndex;
        }
        public static int GetMinElementPosition(double[] array)
        {
            double fMin = array.Min();
            int nIndex = Array.IndexOf(array, fMin);
            return nIndex;
        }

        
    }
}

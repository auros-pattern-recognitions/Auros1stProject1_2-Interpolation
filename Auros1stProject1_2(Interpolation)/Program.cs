using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;

//
// 과제 1-2.
//
// 단위 변환된 파일들의 보간법(Interpolation)
// 파장범위(350~1000)
//
//2021.03.22 조계성


namespace Interpolation
{
    //3중 대각 행렬
    #region TRIDIAGONAL_MATRIX

    // Thomas algorithmn을 이용한 3중 대각 행렬 계산
    // https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
    //
    public sealed class Tridiagonal
    {
        public double[] Solve(double[,] matrix, double[] d)
        {
            int rows = matrix.GetLength(0);
            int cols = matrix.GetLength(1);
            int len = d.Length;
            // 행, 열 개수가 동일한 정방 행렬일 때 동작
            if (rows == cols && rows == len)
            {

                double[] b = new double[rows];//0,0부터의 대각행렬
                double[] a = new double[rows];//1,0부터의 대각행렬
                double[] c = new double[rows];//0,1부터의 대각행렬
                //ax[i-1] + bx[i] + cx[i+1] = d 
                //[b1 c1        ][x1]   = [d1]
                //[a2 b2 c2     ][x2]   = [d2]
                //[   a3 b3 c3  ][x3]   = [d3]
                //[         cn-1][xn-1] = [dn-1]  
                //[        an bn][xn]   = [dn]

                // 3중 대각 행렬의 행과 열 요소 분리
                for (int i = 0; i < rows; i++)
                {
                    for (int j = 0; j < cols; j++)
                    {
                        if (i == j)             //b대각
                            b[i] = matrix[i, j];
                        else if (i == (j - 1))  //c대각
                            c[i] = matrix[i, j];
                        else if (i == (j + 1))  //a대각
                            a[i] = matrix[i, j];
                    }
                }
                try
                {
                    c[0] = c[0] / b[0];
                    d[0] = d[0] / b[0];
                    
                    //앞에서 부터 계산
                    for (int i = 1; i < len - 1; i++)
                    {
                        c[i] = c[i] / (b[i] - a[i] * c[i - 1]);
                        d[i] = (d[i] - a[i] * d[i - 1]) / (b[i] - a[i] * c[i - 1]);
                    }
                    d[len - 1] = (d[len - 1] - a[len - 1] * d[len - 2]) / (b[len - 1] - a[len - 1] * c[len - 2]);

                    // 뒤에서 부터 계산
                    for (int i = (len - 1); i-- > 0;)
                    {
                        d[i] = d[i] - c[i] * d[i + 1];
                    }

                    return d;
                }
                catch (DivideByZeroException)   //b[0]가 0일때 발생하는 예외처리
                {
                    Console.WriteLine("b[0]=0");
                    return null;
                }
            }
            else
            {
                Console.WriteLine("입력한 행렬이 정방행렬이 아닙니다.");
                return null;
            }
        }
    }

    #endregion

    //배열 설정(오름차순, 간격, 스케일
    #region EXTENSIONS
    public static class Extensions
    {
        // 배열을 오름 차순으로 리턴
        public static List<T> SortedList<T>(this T[] array)
        {
            List<T> l = array.ToList();
            l.Sort();
            return l;

        }

        // 배열을 넣으면 배열의 요소의 앞뒤 차이를 배열로 리턴
        public static double[] Diff(this double[] array)
        {
            int len = array.Length - 1;
            double[] diffsArray = new double[len];
            for (int i = 1; i <= len; i++)
            {
                diffsArray[i - 1] = array[i] - array[i - 1];
            }
            return diffsArray;
        }

        // 배열의 값들을 스케일 변환
        public static double[] Scale(this double[] array, double[] scalor, bool mult = true)
        {
            int len = array.Length;
            double[] scaledArray = new double[len];

            if (mult)   //확장
            {
                for (int i = 0; i < len; i++)
                {
                    scaledArray[i] = array[i] * scalor[i];
                }
            }
            else        //축소
            {
                for (int i = 0; i < len; i++)
                {
                    if (scalor[i] != 0)
                    {
                        scaledArray[i] = array[i] / scalor[i];
                    }
                    else
                    {
                        // scalor의 인수가 0일때 예외 발생 방지
                        scalor[i] = 0.00001;
                        scaledArray[i] = array[i] / scalor[i];

                    }
                }
            }
            return scaledArray;
        }
        
        //2차원 행렬에 일차원 배열 대각선으로 저장
        public static double[,] Diag(this double[,] matrix, double[] diagVals)
        {
            int rows = matrix.GetLength(0);// 매트릭스의 크기 참조
            int cols = matrix.GetLength(1);
            // the matrix has to be scare
            if (rows == cols)   //정방행렬
            {
                double[,] diagMatrix = new double[rows, cols];
                int k = 0;
                for (int i = 0; i < rows; i++)
                    for (int j = 0; j < cols; j++)
                    {
                        if (i == j)
                        {
                            diagMatrix[i, j] = diagVals[k]; //일차원 배열을 대각선에 저장
                            k++;
                        }
                        else
                        {
                            diagMatrix[i, j] = 0;
                        }
                    }
                return diagMatrix;
            }
            else
            {
                Console.WriteLine("정방행렬이 아닙니다.");
                return null;
            }


        }
    }

    #endregion

    //3차 곡선 보간법
    #region  CUBIC_SPLINE
    internal sealed class CubicSplineInterpolation
    {
        // f(x) = 1/(1+x^2)*sin(x)

        public CubicSplineInterpolation(double[] _x, double[] _y)
        {
            int xLength = _x.Length;
            if (xLength == _y.Length && xLength > 1)    //배열에 인자 넣기
            {
                x = _x;
                y = _y;
            }
            len = _x.Length;
            if (len > 1)
            {
                baseset = true;
                lX = X.ToList();    //배열을 리스트로 저장
            }
            else
                Console.WriteLine("x,y최소 2개 이상 ");
        }

        public  double? Interpolate(double p)
        {
            if (baseset)
            {
                double? result = 0;
                int N = len - 1;

                double[] h = X.Diff();  //x배열의 차이 저장
                double[] D = Y.Diff().Scale(h, false);//스케일링
                //반복되는 단일 값이 들어있는 시퀀스 생성
                double[] s = Enumerable.Repeat(3.0, D.Length).ToArray();
                //차이값 스케일링
                double[] dD3 = D.Diff().Scale(s);
                double[] a = Y;

                // 3중 대각선 값들 정리
                double[,] H = new double[N - 1, N - 1];
                double[] diagVals = new double[N - 1];
                for (int i = 1; i < N; i++)
                {
                    diagVals[i - 1] = 2 * (h[i - 1] + h[i]);
                }

                H = H.Diag(diagVals);

                // 정사각 행렬이 아닐경우
                if (H != null)
                {
                    for (int i = 0; i < N - 2; i++)
                    {
                        H[i, i + 1] = h[i + 1];
                        H[i + 1, i] = h[i + 1];
                    }

                    double[] c = new double[N + 2];
                    c = Enumerable.Repeat(0.0, N + 1).ToArray();

                    // solve tridiagonal matrix
                    Tridiagonal le = new Tridiagonal();
                    double[] solution = le.Solve(H, dD3);

                    for (int i = 1; i < N; i++)
                    {
                        c[i] = solution[i - 1];
                    }

                    double[] b = new double[N];
                    double[] d = new double[N];
                    for (int i = 0; i < N; i++)
                    {
                        b[i] = D[i] - (h[i] * (c[i + 1] + 2.0 * c[i])) / 3.0;
                        d[i] = (c[i + 1] - c[i]) / (3.0 * h[i]);
                    }

                    double Rx;

                    try
                    {
                        //p값이 x축 좌표계의 범위를 넘어갈때
                        Rx = X.First(m => m >= p);
                    }
                    catch
                    {
                        return null;
                    }

                    //보간 하고자 하는 지점과 가장 가까운 x축 좌표를 찾는다.
                    int iRx = lX.IndexOf(Rx);

                    if (iRx == -1)
                        return null;

                    if (iRx == len - 1 && X[iRx] == p)
                        return Y[len - 1];

                    if (iRx == 0)
                        return Y[0];

                    iRx = lX.IndexOf(Rx) - 1;
                    Rx = p - X[iRx];
                    result = a[iRx] + Rx * (b[iRx] + Rx * (c[iRx] + Rx * d[iRx]));

                    return result;
                }
                else
                {
                    return null;
                }


            }
            else
                return null;
        }
        public double[] X
        {
            get
            {
                return x;
            }
        }

        public double[] Y
        {
            get
            {
                return y;
            }
        }
        private bool baseset = false;
        private int len;
        private List<double> lX;
        private double[] x;
        private double[] y;
    }
    #endregion

    
    class EntryPoint
    {
        static void Interpolation(string path)
        {
           

            string[] MeasurementSpectrumData;   // 측정 스펙트럼 데이터 저장할 배열. (한 줄씩 저장)
            string[] SingleLineData;            // 한 줄의 스펙트럼 데이터를 임시로 저장할 배열.

            MeasurementSpectrumData = File.ReadAllLines(path);
            int LoopNum = MeasurementSpectrumData.Length;

            double[] wavelength = new double[LoopNum - 1];
            double[] n = new double[LoopNum - 1];
            double[] k = new double[LoopNum - 1];

            // 데이터의 첫번째 줄은 column 명이다.
            // 이를 제외하기 위해 반복문을 1부터 시작한다.
            int StartIndex = 1;
            char[] delimiterChars = { ' ', '\t', '\n', ',' };

            for (int i = StartIndex; i < LoopNum; i++)
            {
                // Split한 데이터를 SingleLineData에 저장한다.
                SingleLineData = MeasurementSpectrumData[i].Split(delimiterChars, StringSplitOptions.RemoveEmptyEntries);  

                // 각 컬럼에 해당하는 데이터를 저장한다.
                wavelength[i - 1] = double.Parse(SingleLineData[0]);
                n[i - 1] = double.Parse(SingleLineData[1]);
                k[i - 1] = double.Parse(SingleLineData[2]);
            }
           
            double[] p = new double[131];
            for (int i = 0; i < 131; i++)
            {
                p[i] = 350 + i * 5;
            }

            CubicSplineInterpolation CS1 = new CubicSplineInterpolation(wavelength, n);
            CubicSplineInterpolation CS2 = new CubicSplineInterpolation(wavelength, k);
            
            // 텍스트에 저장
            if (path == @"C:\Users\sksms\Source\Repos\Auros1stProject1_2-Interpolation\Auros1stProject1_2(Interpolation)\data\SiO2_nm.txt")
            {
                path = "SiO2_new.txt";
            }
            else if (path == @"C:\Users\sksms\Source\Repos\Auros1stProject1_2-Interpolation\Auros1stProject1_2(Interpolation)\data\Si_nm.txt")
            {
                path = "Si_new.txt";
            }
            else if (path == @"C:\Users\sksms\Source\Repos\Auros1stProject1_2-Interpolation\Auros1stProject1_2(Interpolation)\data\SiN.txt")
            {
                path = "SiN_new.txt";
            }
            using (StreamWriter NewSpectrumOutputFile = new StreamWriter(path))
            {
                // 컬럼 명 쓰기.
                NewSpectrumOutputFile.WriteLine(
                    $"wavelength(nm)\t" +
                    $"n\t" +
                    $"k");   

                // 스펙트럼 데이터 쓰기.
                foreach (double pp in p)
                {
                    Console.WriteLine($"{pp}\t" +
                        $"{ CS1.Interpolate(pp).ToString() }\t" +
                        $"{ CS2.Interpolate(pp).ToString()}");
                    NewSpectrumOutputFile.WriteLine($"{pp}\t" +
                        $"{ CS1.Interpolate(pp).ToString() }\t" +
                        $" { CS2.Interpolate(pp).ToString()}");
                }
            }


        }
        static void Main(string[] args)
        {
            
            
            string path1 = @"C:\Users\sksms\Source\Repos\Auros1stProject1_2-Interpolation\Auros1stProject1_2(Interpolation)\data\SiO2_nm.txt";
            Interpolation(path1);
            Console.WriteLine();
            string path2 = @"C:\Users\sksms\Source\Repos\Auros1stProject1_2-Interpolation\Auros1stProject1_2(Interpolation)\data\Si_nm.txt";
            Interpolation(path2);
            Console.WriteLine();
            string path3 = @"C:\Users\sksms\Source\Repos\Auros1stProject1_2-Interpolation\Auros1stProject1_2(Interpolation)\data\SiN.txt";
            Interpolation(path3);
           
        }
    }
}
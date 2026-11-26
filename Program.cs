
using System;
using System.Text;

namespace Lab_1_MS
{
    class Program
    {
        static void Main(string[] args)
        {
            Console.OutputEncoding = Encoding.UTF8;

            // --- ВИХІДНІ ДАНІ (ВАРІАНТ 1) ---
            int n = 3; // Кількість вузлів СМО
            int N = 6; // Кількість заявок у замкненій мережі

            // Кількість каналів обслуговування r_i
            int[] r = { 1, 1, 3 };

            // Середні тривалості обслуговування 1 / mu_i
            double[] invMu = { 0.8, 0.3, 1.2 };
            double[] mu = new double[n];
            for (int i = 0; i < n; i++)
                mu[i] = 1.0 / invMu[i];

            // Матриця ймовірностей переходів P (p_13 = 1 - 0.1 - 0.3 = 0.6)
            double[,] P = {
                { 0.1, 0.3, 0.6 },
                { 1.0, 0.0, 0.0 },
                { 1.0, 0.0, 0.0 }
            };

            // 1. Знаходження коефіцієнтів передачі e_i (e_1 = 1)
            double[] e = SolveTransmissionCoefficients(P, n);

            Console.WriteLine("=== КОЕФІЦІЄНТИ ПЕРЕДАЧІ e_i ===");
            for (int i = 0; i < n; i++)
                Console.WriteLine($"e[{i + 1}] = {e[i]:F4}");

            // 2. Обчислення функцій p_i(k) за формулою (1.12)
            double[,] p = new double[n, N + 1];
            for (int i = 0; i < n; i++)
            {
                double rho = e[i] / mu[i];
                for (int k = 0; k <= N; k++)
                {
                    if (k <= r[i])
                    {
                        p[i, k] = Math.Pow(rho, k) / Factorial(k);
                    }
                    else
                    {
                        p[i, k] = Math.Pow(rho, k) / (Factorial(r[i]) * Math.Pow(r[i], k - r[i]));
                    }
                }
            }

            // 3. Обчислення нормуючого множника C(N) та ймовірностей P_CMOi(j)
            double[,] P_CMO = new double[n, N + 1];
            double sumStates = 0.0;

            for (int a = 0; a <= N; a++)
            {
                for (int b = 0; b <= N - a; b++)
                {
                    int c = N - a - b;
                    double stateWeight = p[0, a] * p[1, b] * p[2, c];
                    sumStates += stateWeight;

                    P_CMO[0, a] += stateWeight;
                    P_CMO[1, b] += stateWeight;
                    P_CMO[2, c] += stateWeight;
                }
            }

            double C_N = 1.0 / sumStates;
            Console.WriteLine($"\nНормуючий множник C({N}) = {C_N:F6}\n");

            // Нормуємо ймовірності перебування j заявок у кожній СМО
            for (int i = 0; i < n; i++)
                for (int j = 0; j <= N; j++)
                    P_CMO[i, j] *= C_N;

            // 4. Виведення розподілу ймовірностей P_CMOi(j)
            Console.WriteLine("=== ЙМОВІРНОСТІ СТАНІВ P_CMOi(j) ===");
            Console.Write("j\t");
            for (int i = 0; i < n; i++) Console.Write($"СМО_{i + 1}\t\t");
            Console.WriteLine();

            for (int j = 0; j <= N; j++)
            {
                Console.Write($"{j}\t");
                for (int i = 0; i < n; i++)
                    Console.Write($"{P_CMO[i, j]:F6}\t");
                Console.WriteLine();
            }

            // 5. Обчислення основних характеристик СМО
            Console.WriteLine("\n=== ОСНОВНІ ХАРАКТЕРИСТИКИ МЕРЕЖІ МО ===");
            Console.WriteLine("Вузол\tL_i (черга)\tR_i (зайнято)\tM_i (всього)\tlambda_i\tT_i (перебув.)\tW_i (очікув.)");

            double sumM = 0;
            for (int i = 0; i < n; i++)
            {
                // Середня довжина черги L_i
                double L = 0.0;
                for (int j = r[i] + 1; j <= N; j++)
                    L += (j - r[i]) * P_CMO[i, j];

                // Середня кількість зайнятих каналів R_i
                double idleChannels = 0.0;
                for (int j = 0; j <= r[i] - 1; j++)
                    idleChannels += (r[i] - j) * P_CMO[i, j];
                double R = r[i] - idleChannels;

                // Середня кількість заявок у СМО M_i
                double M = L + R;
                sumM += M;

                // Інтенсивність вихідного потоку lambda_i
                double lambda = R * mu[i];

                // Середній час перебування T_i та чекання W_i
                double T = M / lambda;
                double W = L / lambda;

                Console.WriteLine($"СМО_{i + 1}\t{L:F4}\t\t{R:F4}\t\t{M:F4}\t\t{lambda:F4}\t\t{T:F4}\t\t{W:F4}");
            }

            Console.WriteLine($"\nПеревірка збереження кількості заявок: M_1 + M_2 + M_3 = {sumM:F4} (має дорівнювати N = {N})");
        }

        // Розв'язання системи e_i = sum(e_j * p_ji) при e_1 = 1
        static double[] SolveTransmissionCoefficients(double[,] P, int n)
        {
            // Для n = 3 при e_1 = 1 маємо 2 лінійних рівняння для e_2 та e_3:
            // (1 - p_22)*e_2 - p_32*e_3 = p_12
            // -p_23*e_2 + (1 - p_33)*e_3 = p_13
            double a11 = 1.0 - P[1, 1];
            double a12 = -P[2, 1];
            double b1 = P[0, 1];

            double a21 = -P[1, 2];
            double a22 = 1.0 - P[2, 2];
            double b2 = P[0, 2];

            double det = a11 * a22 - a12 * a21;
            double e2 = (b1 * a22 - b2 * a12) / det;
            double e3 = (a11 * b2 - a21 * b1) / det;

            return new double[] { 1.0, e2, e3 };
        }

        static double Factorial(int k)
        {
            double res = 1.0;
            for (int i = 2; i <= k; i++) res *= i;
            return res;
        }
    }
}
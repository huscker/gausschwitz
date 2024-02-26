/* Голубев Егор ИВТ-32 БО */

#define VERBOSE // закоментировать если нужно отрубить вывод

#include <vector>
#include <iostream>
#include <algorithm>
#include <fstream>

using namespace std;

typedef int datatype;
// Считывание осуществляется из файла `input.txt`
const string filename = "input.txt";




// обмен строк матрицы
void swap_rows(vector<vector<datatype>> &a, int i, int j) {
    vector<datatype> t = a[i];
    a[i] = a[j];
    a[j] = t;
}

// домножение строк на число
void mul_row_by_val(vector<vector<datatype>> &a, int row, datatype val) {
    for (int j = 0; j < a.size(); j++) {
        a[row][j] *= val;
    }
}


// нахождение первой строки у которой ведущий не нулевой
int next_not_null_row(vector<vector<datatype>> &a, int row) {
    for (int i = row + 1; i < a.size(); i++) {
        if (a[i][row] != 0.0) {
            return i;
        }
    }
    return -1;
}


// реализация метода Гаусса
vector<datatype> gauss_elimination(vector<vector<datatype>> A, vector<datatype> b) {
    int n = A.size();

    // Прямой ход
    for (int i = 0; i < n-1; ++i) {
        int not_null_row = next_not_null_row(A, i);
        if (not_null_row == -1) {
            return {}; // Вырожденная система или бесконечное множество решений
        }

        if (not_null_row != i) {
            swap_rows(A, i, not_null_row);
            swap(b[i], b[not_null_row]);
        }

        for (int j = i + 1; j < n; ++j) {
            int gcd = __gcd(abs(A[i][i]), abs(A[j][i]));
            int lcm = (A[i][i] / gcd) * abs(A[j][i]);
            int coeff = lcm / A[j][i];
            mul_row_by_val(A, j, coeff);
            b[j] *= lcm / A[i][i];
            b[j] -= b[i] * coeff;
        }
    }

    // Обратный ход
    vector<datatype> x(n);
    for (int i = n - 1; i >= 0; --i) {
        x[i] = b[i];
        for (int j = i + 1; j < n; ++j) {
            x[i] -= A[i][j] * x[j];
        }
        int gcd = __gcd(abs(A[i][i]), abs(x[i]));
        x[i] /= gcd;
    }

    return x;
}



// Вызов метода Гаусса для считывания матрицы A и вектора b
int main() {
    int n;
    cout << "Считывание осуществляется из файла `input.txt`" << endl;
    ifstream file(filename);
    file >> n;
    vector<vector<datatype>> A(n, vector<datatype>(n));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            file >> A[i][j];
        }
    }
    vector<datatype> b(n);
    for (int i = 0; i < n; ++i) {
        file >> b[i];
    }
    file.close();


    vector<datatype> solution = gauss_elimination(A, b);

    cout << "x:" << endl;
    for (int i = 0; i < n; ++i) {
        cout << solution[i] << " ";
    }

    return 0;
}

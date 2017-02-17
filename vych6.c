#include <stdio.h>
#include <stdlib.h>
#include <math.h>
void copy(double **A, double **arr, int n)
{
    int i, j;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            A[i][j] = arr[i][j];
}
int read(double ** arr)
{
    int n, i, j;
    FILE * f;
    double **x;
    f = fopen("input1.txt", "r");
    fscanf(f, "%d", &n);
    x = malloc(n*sizeof(double*));
    for (i = 0; i < n; i++)
    {
        x[i] = malloc((n + 1)*sizeof(double));
        for (j = 0; j < n + 1; j++)
            fscanf(f, "%lf ", &x[i][j]);
    }
    fclose(f);
    *arr = x;
    return n;
    for (i = 0; i < n; i++)
        free(x[i]);
}
void write(double **arr, int n)
{
    int i, j;
    for (i = 0; i<n; i++)
    {
        printf("\n");
        for (j = 0; j<n; j++)
            printf("%lf ", arr[i][j]);
    }
    printf("\n");
}
double **getunit(int n)
{
    int i;
    double **arr;
    arr = calloc(sizeof(double*), n);
    for (i = 0; i < n + 1; i++)
        arr[i] = calloc(sizeof(double), n);
    for (i = 0; i < n; i++)
        arr[i][i] = 1;
    return arr;
    for (i = 0; i < n; i++)
        free(arr[i]);
}
double ** multimatrix(double **newm, double **old, int n)
{
    int i, j, k;
    double **rezult;
    rezult = getunit(n);
    for (i = 0; i < n; i++)
        rezult[i][i] = 0;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            for (k = 0; k < n; k++)
                rezult[i][j] = rezult[i][j] + newm[i][k] * old[k][j];
    return rezult;
    for (i = 0; i < n; i++)
        free(rezult[i]);
}
double *otrazh(double **arr, int n)
{
    int i, j, m, k;
    double **h = getunit(n);
    double **e = getunit(n);
    double **H = getunit(n);
    double * w = calloc(sizeof(double), n);
    double **aa = getunit(n);
    for (i = 0; i < n - 1; i++)
    {
        h = getunit(n);
        double s = 0;
        for (k = i; k < n; k++)
            s = s + arr[k][i] * arr[k][i];
        double beta = (arr[i][i] < 0) ? sqrt(s) : -sqrt(s);
        double mu = 1 / sqrt(2 * beta*beta - 2 * beta*arr[i][i]);
        w[i] = mu*(arr[i][i] - beta);
        for (j = i + 1; j < n; j++)
            w[j] = mu*arr[j][i];
        for (j = i; j < n; j++)
            for (k = i; k < n; k++)
                h[j][k] = e[j][k] - 2 * w[j] * w[k];
        H = multimatrix(H, h, n);
        for (m = i; m < n; m++)
            for (j = i; j < n + 1; j++)
            {
            aa[m][j] = 0;
            for (k = i; k < n; k++)
                aa[m][j] = aa[m][j] + h[m][k] * arr[k][j];
            }
        for (j = i; j < n; j++)
            for (k = i; k < n + 1; k++)
                arr[j][k] = aa[j][k];
    }

    double s;
    double *c = calloc(sizeof(double), n);
    for (i = n - 1; i >= 0; i--)
    {
        s = arr[i][n];
        for (j = i + 1; j < n; j++)
            s = s - arr[i][j] * c[j];
        c[i] = s / arr[i][i];
    }
    return c;
    for (i = 0; i < n; i++){
        free(h[i]);
        free(e[i]);
        free(aa[i]);
        free(H[i]);
    }
    free(w);
}
double * norm(double *x, int n)
{
    double norm = 0;
    int i;
    for (i = 0; i < n; i++)
        norm = norm + x[i] * x[i];
    norm = sqrt(norm);
    for (i = 0; i < n; i++)
        x[i] = x[i] / norm;
    return x;
}
void print_vec(double *x, int n)
{
    int i;
    for (i = 0; i < n; i++)
        printf("%lf  ", x[i]);
    printf("\n");
}
double lambda(double lam, double *x, double *y, int n)
{
    int i, j,k=0;
    double rezult = 0;
    for (i = 0; i < n; i++)
        if (y[i] != 0)
        {
            rezult = rezult + x[i] / y[i];
            k++;
        }
    rezult = rezult / k + lam;
    return rezult;
}
double **mat(double **B, double lam, double *x, int n)
{
    int i;
    for (i = 0; i < n; i++){
        B[i][i] = B[i][i] - lam;
        B[i][n] = x[i];
    }
    return B;
}
double absol(double x, double y)
{
    if (x - y < 0)
        return -(x - y);
    else return x - y;
}
double obrat(double **arr, int n, double q)
{
    int i, j, k = 0;
    double **B = getunit(n);
    double * x = calloc(sizeof(double), n);
    for (i = 0; i < n; i++)
        x[0] = 1;
    x = norm(x, n);
    double *y = calloc(sizeof(double), n);
    double lam = q, temp = 1000;
    for (i = 0; i < 5; i++){
        copy(B, arr, n);
        B = mat(B, q, x, n);
        y = otrazh(B, n);
        lam = lambda(q, x, y, n);
        x = norm(y, n);
    }
    while (absol(temp, lam)>0.0000001){
        temp = lam;
        copy(B, arr, n);
        B = mat(B, lam, x, n);
        y = otrazh(B, n);
        lam = lambda(lam, x, y, n);
        x = norm(y, n);
        //print_vec(y, n);
        k++;
    }
    print_vec(x, n);
    printf("%d", k);
    return lam;
    for (i = 0; i < n; i++)
        free(B[i]);
    free(x); free(y);
}
int main()
{
    int dim,i;
    double **arr, x, q;
    dim = read(&arr);
    write(arr, dim);
    printf("ENTER PRIBLIZHENIE: ");
    scanf("%lf", &q);
    x = obrat(arr, dim, q);
    printf("\n%lf\n", x);
    return 0;
    for (i = 0; i < dim; i++)
        free(arr[i]);
}



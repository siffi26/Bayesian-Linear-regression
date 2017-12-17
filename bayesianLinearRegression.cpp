/*Program to perform Bayesian Linear regression*/
/*Author: Siffi Singh */
/*Dated: 27/10/2017 */

/* standard Header */
#include <bits/stdc++.h>
using namespace std;
# define two_pi 2.0*3.14159265358979323846
/*Generating guassian distributed data*/
double guassian(double mu, double sigma) {
        double z0 = 0, z1 = 0;
        double epsilon = std::numeric_limits < double > ::min();
        bool generate;
        generate = !generate;
        if (!generate)
            return z1 * sigma + mu;
        double u1, u2;
        do {
            u1 = rand() * (1.0 / RAND_MAX);
            u2 = rand() * (1.0 / RAND_MAX);
        } while (u1 <= epsilon);
        z0 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
        z1 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2);
        return z0 * sigma + mu;
    }
    /*Generating uniform data*/
double uniform(double min, double max) {
        return (min + (rand() % static_cast < int > (max - min + 1)));
    }
    /*matrix multiply for row and column vector*/
double matrixmul(int n, double X[][10], double W[][10]) {
        double ans = 0;
        for (int i = 0; i < 1; i++) {
            for (int j = 0; j < 1; j++) {
                for (int k = 0; k < n; k++) {
                    ans = ans + X[i][k] * W[k][j];
                }
            }
        }
        return ans;
    }
    /*Matrix multiplication*/
void matrixmultiply(int n, int m, int p, int q, double A[10][10], double B[10][10], double res[10][10]) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < q; j++) {
                res[i][j] = 0;
            }
        }
        for (int i = 0; i < n; i++) //no. of rows in the first matrix
        {
            for (int j = 0; j < q; j++) //no. of columns in the second matrix
            {
                for (int k = 0; k < p; k++) //q=m the dimension that is equal
                {
                    res[i][j] += A[i][k] * B[k][j];
                }
            }
        }
    }
    /*Matrix Inverse*/
void inverse(int n, double b[10][10], double a[10][10]) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                a[i][j] = b[i][j];
            }
        }
        for (int i = 0; i < n; i++) {
            for (int j = n; j < 2 * n; j++) {
                if (i == j - n)
                    a[i][j] = 1;
                else
                    a[i][j] = 0;
            }
        }
        for (int i = 0; i < n; i++) {
            double t = a[i][i];
            for (int j = i; j < 2 * n; j++)
                a[i][j] = a[i][j] / t;
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    t = a[j][i];
                    for (int k = 0; k < 2 * n; k++)
                        a[j][k] = a[j][k] - t * a[i][k];
                }
            }
        }
    }
    /*Matrix Transpose*/
void transpose(int n, int m, double a[10][10], double b[10][10]) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                b[i][j] = a[j][i];
            }
        }
    }
    /*Driver Function*/
int main() {
        int n, a;
        /*Generating data as in 1.b*/
        cout << "Enter basis number: ";
        cin >> n;
        cout << "Enter variance of error function: ";
        cin >> a;
        /*Variable Declarations*/
        double min = -10;
        double max = 10;
        double e, x;
        int data;
        n = n + 1;
        double X[1][10], Y[10], WT[10][10], conv_ans[10], E[10], X_data[10][10];
        cout << "Enter the weights: ";
        for (int i = 0; i < n; i++) {
            cin >> WT[i][0]; //Generating transpose of W and creating weight vector
        }
        /*Inputing number of data to be generated.*/
        cout << "Enter the number of data you want to generate: ";
        cin >> data;
        cout << "\n---------------------------------------------\n";
        for (int k = 0; k < data; k++) {
            e = guassian(0, a);	//Guassian error distribution
            x = uniform(min, max);	//Uniform x distribution
            for (int j = n - 1; j >= 0; j--) {
                X[0][j] = pow(x, j);
                X_data[k][j] = X[0][j];
            }
            conv_ans[k] = matrixmul(n, X, WT);
            E[k] = e;
            Y[k] = matrixmul(n, X, WT) + e;
            cout << Y[k] << " || " << x << " || " << e << endl;
        }
        double var_y, mean_y, b;
        cout << "\nEnter the precision: ";
        cin >> b;
        b = 1 / b; //Precision = inverse of variance
        var_y = a; //variance of predictive distribution y
        double xct[10][10], xc[10][10], xtx[10][10] = {
            0
        };
        /*calculating lambda or variance of the posterior distribution*/
        for (int k = 0; k < data; k++) {
            for (int j = 0; j < n; j++) {
                xc[0][j] = X_data[k][j];	//Creating vector x
                xct[j][0] = X_data[k][j];	//creating vector X(transpose)
            }
            matrixmultiply(n, 1, 1, n, xct, xc, xtx);

            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    xtx[i][j] *= a;	//Computing X*X(transpose) for lambda
                }
            }
            double bI[10][10] = {
                0
            };
            /*Creating b*I identity matrix*/
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    if (i == j) {
                        bI[i][j] = b;
                    }
                }
            }
            /*Adding a*XT*X and bI matrix to form a*XT*X + bI*/
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    xtx[i][j] += bI[i][j];
                }
            }
            cout<<"\nNew Mean and New variance are:"<<endl;
			/*printing variance of posterior distribution which is the matrix xtx[][]*/
			cout<<"\nvariance of posterior distribution:"<<endl;
			for(int i=0; i<n; i++)
			{
				for(int j=0; j<n; j++)
				{
					cout<<xtx[i][j]<<" "; 
				}
				cout<<endl;
			}
            /*lambda, i.e the variance of the mutivariate guassian w, is xtx[][] matrix*/
            /*finding inverse of lambda matrix to compute mean of the mutivariate guassian*/
            /*Mean of posterior(mu) : a*lambda(inverse)*X(transpose)*y and Variance of posterior(lambda): a*X(transpose)*X + b*I*/
            double lambda_inv[10][10];
            inverse(n, xtx, lambda_inv); //Variance computed
            /*Computing mean of posterior distribution*/
            double mean_pos[10][10];
            matrixmultiply(n, n, n, 1, lambda_inv, xct, mean_pos);
            double mean_sum = 0, w_mean = 0;
            /*printing variance of posterior distribution which is the vector mean_pos[]*/
            cout<<"\nmean of posterior distribution:"<<endl;
            for (int i = 0; i < n; i++) {
            	cout<<mean_pos[i][0]<<" ";
                mean_sum += mean_pos[i][0];
                w_mean += WT[i][0];
            }
            cout<<endl;
            mean_sum /= n; //mean of mean vector of posterior distribution
            w_mean /= n; //mean of w
            //cout<<mean_sum<<" "<<w_mean<<endl;

            /*calculating mean and variance of univariate guassian of Predictive Distribution*/
            /*Mean of predictive : mu(transpose)*X(transpose) and Variance of predictive: (1/a)*X(transpose)*lambda(inverse)*X. */

            transpose(1, n, mean_pos, mean_pos); //transpose of mean for mu(T), mean of predictive distribution
            double mean_pred = matrixmul(n, mean_pos, xct);

            cout<<"\nmean of predictive distribution:"<<endl;
            cout << mean_pred << endl; //Printing mean of the predictive distribution 'y'
            double var_pred[10][10], ans[10][10], var_p;

            matrixmultiply(1, n, n, n, xc, lambda_inv, var_pred);
            matrixmultiply(1, n, n, 1, var_pred, xct, ans);
            var_p = ans[0][0] + (1 / a);

            cout<<"\nvariance of predictive distribution:"<<endl;
            cout << ans[0][0] + (1 / a) << endl; //Printing variance of the predictive distribution 'y'

            if (mean_sum == w_mean) //Convergence of mean
                cout<<"Mean is converged after "<<k<<"th iteration"<<endl;
            if(var_y == var_p)	//Convergence of variance
                cout<<"Variance is converged after "<<k<<"th iteration"<<endl;
        }
	
	return 0;
}




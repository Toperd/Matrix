#include <iostream>

using namespace std;

//paste.ubuntu.com/15333449/


class Matrix
{
protected:
	int n; // ñòðîêè
	int m; // ñòîëáöû
	double* data;
public:
	Matrix()
	{
		this->n = 0;
		this->m = 0;
		this->data = NULL;
	}
	Matrix(int m, int n, double* data)
	{
		this->n = n;
		this->m = m;
		this->data = data;
	}
	~Matrix()
	{
		delete[] this->data;
	}
	void set(int j, int i, double data)
	{
		if(i > this->n - 1)
		{
			i = this->n - 1;
		}
		if(i < 0)
		{
			i = 0;
		}
		if(j > this->m - 1)
		{
			j = this->m - 1;
		}
		if(j < 0)
		{
			j = 0;
		}
		this->data[(i * this->m) + j] = data;
	}
	double get(int j, int i)const
	{
		return this->data[i * this->m + j];
	}
	int getM()
	{
		return this->m;
	}
	int getN()
	{
		return this->n;
	}
	Matrix(const Matrix& A)
	{
		this->m = A.m;
		this->n = A.n;
		if(A.data != NULL)
		{
			this->data = new double[this->n * this->m];
			for(int i = 0; i < this->n; i++)
			{
				for(int j = 0; j < this->m; j++)
				{					
					this->set(j, i, A.get(j, i));
				}
			}
		}
		else
		{
			this->data = NULL;
		}
	}
	bool failed()
	{
		return(this->data == NULL || this->m == 0 || this->n == 0);
	}
	Matrix& operator = (const Matrix& A)
	{
		if(this->data != NULL)
		{
			delete[] this->data;
		}
		this->m = A.m;
		this->n = A.n;
		this->data = new double[this->n * this->m];
		for(int i = 0; i < this->n; i++)
		{
			for(int j = 0; j < this->m; j++)
			{
				this->set(j, i, A.get(j, i));
			}
		}
		return *this;
	}
	Matrix operator + (Matrix& A)
	{
		Matrix result;
		result.n = 0;
		result.m = 0;
		result.data = NULL;
		if(this->n == A.n && this->m == A.m && !(this->failed()) && !(A.failed()))
		{
			result.n = A.n;
			result.m = A.m;
			result.data = new double [result.n * result.m];
			for(int i = 0; i < A.n; i++)
			{
				for(int j = 0; j < A.m; j++)
				{
					result.set(j, i, (this->get(j, i) + A.get(j, i)));
				}
			}			
		}
		return Matrix(result);
	}
	Matrix operator * (Matrix& A)
	{
		Matrix result;
		result.n = 0;
		result.m = 0;
		result.data = NULL;
		if(!(this->failed()) && !(A.failed()) && this->m == A.n)
		{
			result.n = this->n;
			result.m = A.m;
			result.data = new double [result.n * result.m];
			for(int i = 0; i < this->n; i++)
			{
				for(int j = 0; j < A.m; j++)
				{
					double tmp = 0;
					for(int k = 0; k < A.n; k++)
					{
						tmp += this->get(k, i) * A.get(j, k);
					}
					result.set(j, i, tmp);
				}
			}
		}
		return Matrix(result);
	}
	Matrix operator * (double& C)
	{
		Matrix result;
		result.n = 0;
		result.m = 0;
		result.data = NULL;
		if(!(this->failed()))
		{
			result.n = this->n;
			result.m = this->m;
			result.data = new double [result.n * result.m];
			for(int i = 0; i < result.n; i++)
			{
				for(int j = 0; j < result.m; j++)
				{
					result.set(j, i, C * this->get(j, i));
				}
			}
			return Matrix(result);
		}
	}
	Matrix operator - (Matrix& A)
	{
		Matrix result;
		result.n = 0;
		result.m = 0;
		result.data = NULL;
		if(this->n == A.n && this->m == A.m && !(this->failed()) && !(A.failed()))
		{
			result.n = A.n;
			result.m = A.m;
			result.data = new double [result.n * result.m];
			for(int i = 0; i < A.n; i++)
			{
				for(int j = 0; j < A.m; j++)
				{
					result.set(j, i, (this->get(j, i) - A.get(j, i)));
				}
			}			
		}
		return Matrix(result);
	}
	Matrix transpose()
	{
		Matrix result;
		result.n = 0;
		result.m = 0;
		result.data = NULL;
		if(!(this->failed()))
		{
			result.n = this->m;
			result.m = this->n;
			result.data = new double[result.n * result.m];
			for(int i = 0; i < this->n; i++)
			{
				for(int j = 0; j < this->m; j++)
				{
					result.set(i, j, this->get(j, i));
				}
			}
		}
		return Matrix(result);
	}
	Matrix minor_(int j, int i)
	{
		Matrix result;
		result.n = 0;
		result.m = 0;
		result.data = NULL;
		if(this->n != 1 && this->m != 1 && !(this->failed()))
		{
			result.n = this->n - 1;
			result.m = this->m - 1;
			result.data = new double[result.n * result.m];
			for(int k = 0; k < i; k++)
			{
				for(int p = 0; p < j; p++)
				{
					result.set(p, k, this->get(p, k));
				}
				for(int s = j + 1; s < this->m; s++)
				{
					result.set(s - 1, k, get(s, k));
				}
			}
			for(int u = i + 1; u < this->n; u++)
			{
				for(int q = 0; q < j; q++)
				{
					result.set(q, u - 1, this->get(q, u));
				}
				for(int t = j + 1; t < this->m; t++)
				{
					result.set(t - 1, u - 1, get(t, u));
				}
			}
		}
		return Matrix(result);
	}
	double determinant()
	{
		double result = 0;
		if(this->n == this->m && !(this->failed()))
		{
			if(this->n == 1)
			{
				return this->get(0, 0);
			}
			else
			{				
				int i = 0;
				for(int j = 0; j < this->m; j++)
				{
					if((i + j) % 2 == 0)
					{
						result = result + this->get(j, i) * (this->minor_(j, i)).determinant();
					}
					else
					{
						result = result - this->get(j, i) * (this->minor_(j, i)).determinant();
					}					
				}
				return result;
			}
		}
		else
		{
			return 350589;
		}
	}
	Matrix reverse()
	{
		Matrix result;
		result.n = 0;
		result.m = 0;
		result.data = NULL;
		if(this->determinant() != 0 && this->m == this->n && !(this->failed()))
		{
			result.n = this->n;
			result.m = this->m;
			result.data = new double[result.n * result.m];
			double det = 1 / ( 1.0 * this->determinant());
			for(int i = 0; i < this->n; i++)
			{
				for(int j = 0; j < this->m; j++)
				{
					if((i + j) % 2 == 0)
					{
						result.set(j, i, this->minor_(j, i).determinant());
					}
					else
					{
						result.set(j, i, (-1) * (this->minor_(j, i).determinant()));
					}
				}
			}
			result = result.transpose();
			result = result * det;
		}
		return Matrix(result);

	}
	ostream& print(ostream& o)
	{
		for(int i = 0; i < this->n; i++)
		{
			for(int j = 0; j < this->m; j++)
			{
				o << this->get(j, i) << '\t';
			}
			o << '\n';
		}
		return o;
	}

	istream& read(istream& in)
	{
		double tmp;
		in >> this->m;
		in >> this->n;
		this->data = new double[this->m * this->n];
		for(int i = 0; i < this->n; i++)
		{
			for(int j = 0; j < this->m; j++)
			{
				in >> tmp;
				this->set(j, i, tmp);
			}
		}
		return in;
	}
};

Matrix* get_init(int m, int n) // âîçâðàùàåò óêàçàòåëü íà ìàòðèöó m íà n ñ íóëåâûì çàïîëíåíèåì
{
	double* mass = new double[n * m];
	for(int i = 0; i < n * m; i++)
	{
		mass[i] = 0;
	}
	Matrix* test = new Matrix(m, n, mass);
	return test;
}
/*int main()
{

	Matrix* test = get_init(3, 3);

	test->print(cout);

	test->set(0, 0, 2);
	test->set(0, 1, 3);
	test->set(0, 2, 4);
	test->set(1, 0, 1);
	test->set(1, 1, 1);
	test->set(1, 2, 0);
	test->set(2, 0, 0);	
	test->set(2, 1, 0);
	test->set(2, 2, 1);
	

	cout << '\n';
	test->print(cout);
	
	Matrix test2;

	cout << '\n' << test->determinant() << '\n';

	test2 = test->transpose();
	test2.print(cout);

	if(test->failed())
	{
		cout << "failed" << '\n';
	}

	return 0;
}*/

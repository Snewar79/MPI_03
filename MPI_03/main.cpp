


#include <mpi.h>
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <vector>
#define N 3
#define min_size 3


using namespace std;


struct  point;
void ConvexHullJarvis(const vector<point> &mas, vector<int> &convex_hull, int n);
bool More(double x, double y);
bool Less(double x, double y);
bool Equal(double x, double y);
double dist(point, point);
double CosAngle(point a, point b, point c);

struct point
{
	int x;
	int y;

	point(int _x, int _y)
	{
		x = _x;
		y = _y;
	}
	point()
	{
	
	};

};


bool operator!=(point a, point b)
{
	if (a.x != b.x || a.y != b.y) return true;
	return false;
}

bool operator==(point a, point b)
{
	if (a.x == b.x && a.y == b.y) return true;
	return false;
}






double v_cos(int x0, int y0, int x1, int y1, int x2, int y2)
{
	int ax, ay, bx, by;

	ax = x1 - x0;
	ay = y1 - y0;
	bx = x2 - x1;
	by = y2 - y1;

	double cos;

	cos = ((ax * bx) + (ay * by)) / (sqrt(ax * ax + ay * ay) * sqrt(bx * bx + by * by));
	return cos;
}



void shell(int *buffer, int point_count, int **ret_buf)
{
	point p0;
	point p1;
	point p2;
	vector<point> v_p;
	vector<int> out_index;

	point *point_buffer = new point[point_count];
	point *shell = new point[point_count];

	p1.x = buffer[0];
	p1.y = buffer[1];

	int shell_count = 0;

	for (int i = 0; i < point_count; i++)
	{
		if (buffer[i * 2 + 1] < p1.y)
		{
			p1.x = buffer[i * 2];
			p1.y = buffer[i * 2 + 1];
		}
		point_buffer[i].x = buffer[i * 2];
		point_buffer[i].y = buffer[i * 2 + 1];
		v_p.push_back(point_buffer[i]);
	}

	ConvexHullJarvis(v_p, out_index, point_count);

	for (int i = 0; i < out_index.size(); i++)
	{
		cout << out_index[i] << "\n";
	}
	
}





void ConvexHullJarvis(const vector<point> &mas, vector<int> &convex_hull, int n)
{
	// находим самую левую из самых нижних
	int base = 0;
	for (int i = 1; i<n; i++)
	{
		if (mas[i].y < mas[base].y)
			base = i;
		else
			if (mas[i].y == mas[base].y &&
				mas[i].x < mas[base].x)
				base = i;
	}
	// эта точка точно входит в выпуклую оболочку
	convex_hull.push_back(base);

	point first = mas[base];
	point cur = first;
	point prev = point(first.x - 1, first.y);
	do
	{
		double minCosAngle = 1e9; // чем больше угол, тем меньше его косинус
		double maxLen = 1e9;
		int next = -1;


		for (int i = 0; i<n; i++)
		{
			
			double curCosAngle = CosAngle(prev, cur, mas[i]);
			if (Less(curCosAngle, minCosAngle))
			{
				next = i;
				minCosAngle = curCosAngle;
				maxLen = dist(cur, mas[i]);
			}
			else if (Equal(curCosAngle, minCosAngle))
			{
				double curLen = dist(cur, mas[i]);
				if (More(curLen, maxLen))
				{
					next = i;
					maxLen = curLen;
				}
			}
		}



		prev = cur;
		cur = mas[next];
		convex_hull.push_back(next);
	} while (cur != first);
}


bool More(double x, double y)
{
	if (x > y) return true;
	return false;
}

bool Less(double x, double y)
{
	if (x < y) return true;
	else return false;
}


bool Equal(double x, double y)
{
	if (x == y) return true;
	else return false;
}


double dist(point a, point b)
{
	return sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
}


double CosAngle(point a, point b, point c)
{

	if (a == b || b == c || c == a) return 1e10;

	int ax, ay, bx, by;

	ax = b.x - a.x;
	ay = b.y - a.y;
	bx = c.x - b.x;
	by = c.y - b.y;

	double cos;

	cos = ((ax * bx) + (ay * by)) / (sqrt(ax * ax + ay * ay) * sqrt(bx * bx + by * by));
	return cos;
	

}


int main(int argc, char **argv)
{

	int taskid, numtasks;

	int work_size;

	int *point_buffer = NULL;

	int *p_work_size = NULL;

	int *disp = NULL;

	int *double_p_work_size = NULL;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);


	if (taskid == 0)
	{

		srand(time(NULL));

		work_size = N;
		
		point_buffer = new int[2 * N];

		/*for (int i = 0; i < N * 2; i++)
		{
			point_buffer[i] = rand() % 100 - 50;
		}
		*/



		point_buffer[0] = 3;
		point_buffer[1] = 0;
		point_buffer[2] = 4;
		point_buffer[3] = 3;
		point_buffer[4] = -3;
		point_buffer[5] = 2;





		cout << "Process NUM = " << numtasks << "\n";

		string points = "";

		for (int i = 0; i < N; i++)
		{
			points += "(" + to_string(point_buffer[i * 2]) + ", " + to_string(point_buffer[i * 2 + 1]) + ") ";
		}
		cout << points << "\n";



		p_work_size = new int[numtasks];
		
		int standart_task;

		standart_task = N / numtasks;

		cout << "standart_task = " << standart_task << "\n";
		
		int real_proc_count = 0;

		if (standart_task >= min_size) // All is ok. «аданий достаточно
		{
			for (int i = 0; i < numtasks - 1; i++)
			{
				p_work_size[i] = standart_task;
			}
			p_work_size[numtasks - 1] = N % numtasks + standart_task;

			real_proc_count = numtasks;
		}
		else // “очек недостаточно
		{
			real_proc_count = N / min_size;
			if (N % min_size != 0) real_proc_count++;

			for (int i = 0; i < real_proc_count; i++)
			{	
				p_work_size[i] = min_size;
			}
			
			if (N % min_size != 0)
			p_work_size[real_proc_count - 1] = N % min_size;
			else
			{
				p_work_size[real_proc_count - 1] = min_size;
			}

			for (int i = real_proc_count; i < numtasks; i++)
			{
				p_work_size[i] = 0;
			}

		}

		// ѕосчитали, значит количество точек и минимальные задани€
		for (int i = 0; i < numtasks; i++)
		{
			cout << i << " = " << p_work_size[i] << " | ";
		}
		cout << "\n";
		// «ашибись все поделено. ќсталось разослать
		//~
		// јн нет, надо бы еще смещени€ посчитать

		disp = new int[numtasks];
		disp[0] = 0;

		for (int i = 1; i < numtasks; i++)
		{
			disp[i] = disp[i - 1] + p_work_size[i - 1] * 2;
		}

		double_p_work_size = new int[numtasks];


		for (int i = 0; i < numtasks; i++)
		{
			double_p_work_size[i] = p_work_size[i] * 2;
		}

	}

	MPI_Scatter(p_work_size, 1,	MPI_INT, &work_size, 1,	MPI_INT, 0,	MPI_COMM_WORLD);

	//cout << "task id = " << taskid << " work_size " << work_size << "\n";

	int *rec_point_buf = new int[work_size * 2];

	MPI_Scatterv(point_buffer, double_p_work_size, disp, MPI_INT, rec_point_buf, work_size * 2, MPI_INT, 0, MPI_COMM_WORLD);

	/*string points = to_string(taskid) + "points ";
	for (int i = 0; i < work_size; i++)
	{
		points += "(" + to_string(rec_point_buf[i * 2]) + ", " + to_string(rec_point_buf[i * 2 + 1]) + ") ";
	}
	cout << points << "\n";
	*/
	int *shell_buf;

	shell(rec_point_buf, work_size, &shell_buf);

	

	MPI_Finalize();

	return 0;
}



#include <mpi.h>
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <vector>
#define N 10
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







int shell(int *buffer, int point_count, vector<int> &out_index)
{
	
	vector<point> v_p;
	
	int i_p1 = 0;
	int i_p2 = -1;
	int i_p3 = -1;
	// ѕосчитать количество различных точек. »зи

	
	point *shell = new point[point_count];

	for (int i = 0; i < point_count; i++)
	{
		point point_buffer;
		point_buffer.x = buffer[i * 2];
		point_buffer.y = buffer[i * 2 + 1];
		v_p.push_back(point_buffer);

		if (v_p[i] != v_p[i_p1])
		{
			if (i_p2 == -1)
			{
				i_p2 = i;
			}
			else
			{
				if (v_p[i_p2] != v_p[i])
				{
					i_p3 = i;
				}

			}

		}
	}

	if (i_p1 == 0 && i_p2 == -1 && i_p2 == -1)
	{
		out_index.push_back(0);
		return 1;
	}
	if (i_p1 == 0 && i_p2 != -1 && i_p3 == -1)
	{
		out_index.push_back(0);
		out_index.push_back(i_p2);
		return 2;
	}


	ConvexHullJarvis(v_p, out_index, point_count);


	
	int shell_count = out_index.size() - 1;
	

	return shell_count;

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
		double minCosAngle = -1e9; // чем больше угол, тем меньше его косинус
		double maxLen = 1e9;
		int next = -1;


		for (int i = 0; i<n; i++)
		{
			
			double curCosAngle = CosAngle(prev, cur, mas[i]);
			if (More(curCosAngle, minCosAngle))
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

	if (a == b || b == c || c == a) return -1e10;

	int ax, ay, bx, by;

	ax = b.x - a.x;
	ay = b.y - a.y;
	bx = c.x - b.x;
	by = c.y - b.y;

	a.x *= -1;
	a.y *= -1;

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

		for (int i = 0; i < N * 2; i++)
		{
			//point_buffer[i] = rand() % 100 - 50;
			point_buffer[i] = 3;
		}
		
		point_buffer[N * 2 - 1] = -3;

	



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
	int *shell_buf = NULL;

	int shell_count;

	vector<int> out_index;

	shell_count = shell(rec_point_buf, work_size, out_index);


	shell_buf = new int[shell_count * 2];

	for (int i = 0; i < shell_count; i++)
	{
		shell_buf[i * 2] = rec_point_buf[out_index[i] * 2];
		shell_buf[i * 2 + 1] = rec_point_buf[out_index[i] * 2 + 1];
	}

	/*string points = to_string(taskid) + "-id points " + "diff = " + to_string(work_size - shell_count) + " ";

	for (int i = 0; i < shell_count; i++)
	{
	points += "(" + to_string(shell_buf[i * 2]) + ", " + to_string(shell_buf[i * 2 + 1]) + ") ";
	}
	cout << points << "\n"; */
	
	int *shell_count_buf = NULL;

	if (taskid == 0) shell_count_buf = new int[numtasks];


	MPI_Gather(&shell_count, 1, MPI_INT, shell_count_buf, numtasks, MPI_INT, 0, MPI_COMM_WORLD);

	int * last_shell_buf;
	int *recvcounts;
	int last_shell_count;

	if (taskid == 0)
	{
		last_shell_count = 0;
		for (int i = 0; i < numtasks; i++)
		{



		}
	}

	MPI_Gatherv(shell_buf, shell_count * 2, MPI_INT,
		last_shell_buf, recvcounts, disp, MPI_INT,
		0, MPI_COMM_WORLD);

	MPI_Finalize();

	return 0;
}
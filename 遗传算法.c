#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctime>
#include <math.h>

#define N 30000000000
#define PI 3.14159265
#define MIN_(a,b) ((a)<(b)?(b):(a))//定义求最大值函数MIN_(a,b) 

#define SIZE  50
#define MAXGEN  50
#define P_CORSS 0.75
#define P_MUTATION 0.05

#define LEN 22

typedef struct node//定义结构体，表示染色体，x是染色体二进制编码，fitness是该染色体的适应值，fitsum是当前所有染色体的适应值的总和
{
	char x[LEN];
	double fitness, fitsum;
}node;

node cur[SIZE], next[SIZE], max, min;//定义结构体数组，有当前染色体种群，下一代染色体种群，以及最优和最差染色体。

double randd()
{
	return (double)rand() / RAND_MAX;//产生0-1之间随机数， RAND_MAX其值至少为32767，rand()在0-MAX之间
}
int randi(int k)
{
	return (int)(randd()*k);//返回接近0-k之间的随机值
}

//计算当前种群中各个个体的适应度 
void cal_fitness()
{
	int i, j, k;
	double d;
	for (i = 0; i<SIZE; i++)
	{
		k = 0;
		for (j = LEN - 1; j >= 0; j--) k = (k << 1) + cur[i].x[j];;//该循环用于将染色体二进制数转成十进制，<<是左移运算，x[j]是染色体i的第j位
		d = (double)k / N;//用上面的k产生一个浮点数d，作为函数自变量x的值
		cur[i].fitness = 7 * cos(4* d);            //计算每个染色体的适应值  函数y= 7*cos(4*x)//
		cur[i].fitsum = i>0 ? (cur[i].fitness + cur[i - 1].fitsum) : (cur[0].fitness);//计算所有适应值的总和，里面有迭代过程，就是若i=3,则cur[3].fitsum=cur[3].fitsum+cur[2].fitsum+cur[1].fitsum+cur[0].fitsum.
	}
}

void init()//初始化染色体种群，共SIZE个，每个都是二进制LEN长度的串
{
	int tmp;
	for (int i = 0; i<SIZE; i++)
	{
		tmp = randi(N);//用randi函数产生一个大概N那么大的随机数
		for (int j = 0; j<LEN; j++)
		{
			cur[i].x[j] = tmp % 2;//用temp对2取余，初始化每个染色体的每个位，共LEN位;一个数不是奇数就是偶数，对2取余得0或1，从而产生像001001010...11010001这样的串
			tmp = tmp >> 1;//右移tmp，就是将temp缩小，使tmp值变化
		}
	}
	cal_fitness();//计算当前种群各个个体的适应度 
}

int sel()//选择种群中某个染色体
{
	double p = randd();
	double sum = cur[SIZE - 1].fitsum;//全部染色体的适应值总和
	for (int i = 0; i<SIZE; i++)
	{
		if (cur[i].fitsum / sum>p) return i;//如果当前适应总值/全部总值大于随机数p，也就是当前适应总值不为0，不太差时，返回这个i值，循环停止
	}
}

//换代 
void tran()
{
	int i, j, pos;
	//找当前种群最优个体 
	max = cur[0];
	for (i = 1; i<SIZE - 1; i++)
	{
		if (cur[i].fitness>max.fitness)  max = cur[i];//记录种群中的最优染色体为max
	}
	for (int k = 0; k<SIZE; k += 2)
	{
		//选择交叉个体 
		i = sel();//用sel函数挑选一个个体，染色体号为i
		j = sel();//同上，再选择个个体，号为j

		//选择交叉位置 
		pos = randi(LEN - 1);//随机产生个LEN-1左右的数pos，如18,19,21等

		//交叉
		if (randd()<P_CORSS)//如果randd()产生的随机数(在0-1之间）小于设定的交叉率P_CORSS
		{
			memcpy(next[k].x, cur[i].x, pos);//提取染色体i的前pos位赋给下一代种群next的第k个染色体.——函数memcpy（&a,&b,n)用于从&b的位置开始数n个长度的数据拷贝赋给&a; 
			memcpy(next[k].x + pos, cur[j].x + pos, LEN - pos);//提取染色体j的 后面LEN-pos位赋给next的第k个染色体,结合上面，从而拼成一个新的染色体

			memcpy(next[k + 1].x, cur[j].x, pos);//同样的方式给next的第k+1染色体赋值，这回换过来生成这个值，即提取j的前pos位数据 + i的后LEN-pos位数据
			memcpy(next[k + 1].x + pos, cur[i].x + pos, LEN - pos);
		}
		else//否则不交叉，
		{
			memcpy(next[k].x, cur[i].x, LEN);
			memcpy(next[k + 1].x, cur[j].x, LEN);
		}
		//变异
		if (randd()<P_MUTATION)//如果一随机数小于设定的变异率P_MUTATION，则执行变异操作
		{
			pos = randi(LEN - 1);//仍然是从中间找个位置，执行变异
			next[k].x[pos] ^= next[k].x[pos];// ^=按位异或后赋值函数，相同为0不同为1，此处是将第pos位上的值（不论是0或1）定为0

			pos = randi(LEN - 1);
			next[k + 1].x[pos] ^= next[k + 1].x[pos];//同样将下一代next的第k+1个染色体的第pos位也变为0；？？最好应该是0变1,1变0...,但这样也行，进化稍慢点
		}
	}
	//找下一代的最差个体 
	min = next[0], j = 0;
	for (i = 1; i<SIZE - 1; i++)
	{
		if (next[i].fitness<min.fitness)  min = next[i], j = i;//用j记录最差染色体号
	}
	//用上一代的最优个体替换下一代的最差个体
	next[j] = max;

	memcpy(cur, next, sizeof(cur));//把整个改良过的（经过交叉，变异，替换等操作的）next代的值赋给当前种群cur，供下一次循环优化


	cal_fitness();
}

//打印个体适应度和二进制编码 
void print(node tmp)
{
	printf("%.6lf", tmp.fitness);
	for (int i = 0; i<LEN; i++)  printf(" %d", tmp.x[i]);
	printf("\n");
}

//打印种群
void printcur()
{
	for (int i = 0; i<SIZE; i++) print(cur[i]);
}


void GA()
{
	int cnt = 0;
	double ans;
	while (cnt++<MAXGEN)//当计数值cnt小于设定的最大进化次数MAXGEN，执行换代操作
	{
		tran();

		//    printf("%.6lf\n",max.fitness);
		//    printcur();
	}
	ans = cur[0].fitness;
	for (int i = 1; i<SIZE; i++) ans = MIN_(ans, cur[i].fitness);//找出函数最小值，打印输出（应设定的是求函数最小值）
	printf("%.6lf\n", ans);
}

int main()
{
	srand((unsigned)time(NULL));//初始化随机数种子产生器，为使程序中每次的rand()函数产生的随机数不一样

	init();//初始化种群
	GA();//遗传换代操作

	system("pause"); //输出结果在屏幕，而不是一闪而过
	return 0;
}

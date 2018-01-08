//-------------------------------------------------by Ye Liang at 18:38 on 16th Dec. 2014------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define NMAX 800
#define MMAX (NMAX+31)/32
#define MAX_LOOP 5 //this is just a coefficient of the maximum loop limitation: maximum loop number = MAX_LOOP (times) # of vertices(Scale of the problem)[Adaptive Design]
#define GROUPSIZE 7//this is the groupsize, after selection only this number of individuals will be there. and only these individual could be able to reproduce new ones.
#define CROSSOVERSIZE (int)((((float)GROUPSIZE + 1) / 2)*(GROUPSIZE ))//means the maximum possible individuals after a reproduction process(before selection).
#define FENCESIZE 100 // the size of fences(I use in the crossover process)
#define SETWD(pos) ((pos)>>5)//define basic computation operation
#define SETBT(pos) ((pos)&037)//define basic computation operation
#define ADD_ELEMENT(setadd, pos) ((setadd)[SETWD(pos)] |= bit[SETBT(pos)])//define basic computation operation: add 1 to assigned position
#define DEL_ELEMENT(setadd, pos) ((setadd)[SETWD(pos)] &= ~bit[SETBT(pos)])//define basic computation operation: minus 1 to assigned position
#define IS_ELEMENT(setadd, pos) ((setadd)[SETWD(pos)] &  bit[SETBT(pos)])//define basic computation operation: examine if it is 1 in the assigned position.
#define POP_COUNT(x) (bytecount[(x) >> 24 & 0377] + bytecount[(x) >> 16 & 0377] + bytecount[(x) >> 8 & 0377] + bytecount[(x) & 0377])
int bit[]=//define bit used for position calculation.
{	       020000000000, 010000000000,
04000000000, 02000000000, 01000000000,
 0400000000,  0200000000,  0100000000,
  040000000,   020000000,   010000000,
   04000000,    02000000,    01000000,
    0400000,     0200000,     0100000,
	 040000,      020000,      010000,
	  04000,       02000,       01000,
	   0400,        0200,        0100,
	    040,         020,         010,
		 04,          02,          01
};
int bytecount[]=
{0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,
1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8};
int nCprEdg;
int nVerNum;
int random[1000000];
int Fence[FENCESIZE];//this is a group of 32-bit 0-1 fences used in crossover process, 1 means this part is from the father, and 0 means this part is from the mother, 
//the # of 0s and 1s are equal in statistics. These fences can be used to simulate random crossover.
int fence_i = 0;//fence index.
int random_i = 0;//random number index.
//I always like to generate these fences and random numbers before the main algorithms runs, and whenever I need them, I just pick one. I believe it saves time.

int randf(void);//to pick up a previously generated random integer.

void makerand(void);//to make 1000000 random integers. I do this at the begining of the program, and whenever I need an random integer, I just pick one.

int Set_Size(int Set[]);//to count the number of ones for a compressed adjacency matrix data set

int Read_Graph(int G[NMAX][MMAX]);//to read the graph information into a compressed adjacency matrix.

int Print_Set(int n, int Set[], int j);//to print the data set for a compressed adjacency matrix data set.

int Can_Be_Blue(int G[NMAX][MMAX], int Blue[]);//to test if this point can be excepted from the dominating set.

int Cross_Over(int Generation[CROSSOVERSIZE][NMAX]);//The cross over process in genetic algorithm

int Mutation(int G[NMAX][MMAX], int Generation[CROSSOVERSIZE][NMAX]);//the mutation process in genetic algorithm

int BFS_Array(int BfsArray[NMAX], int G[NMAX][MMAX], int Red[MMAX], int Root);//Using BFS search to generate a SearchOrder

int BinaryG_to_Generaiton(int G[NMAX][MMAX], int BinaryGen[], int Generation[]);//Translate the domset from Compressed Adj Matrix to Genetic Chromosome Data Structure

int Survival_Of_The_Fittest(int G[NMAX][MMAX], int Generation[CROSSOVERSIZE][NMAX]);//Evaluation, ranking, and selection process is in this function.

int Initialise(int G[NMAX][MMAX], int Generation[CROSSOVERSIZE][NMAX], int DomtedSet[]);//initialise the individuals with the principles in the exact algorithm

void Generate_Dom_Set(int G[NMAX][MMAX], int Red[], int Blue[], int DomtingSet[], int DomtedSet[], int nDepthIndex,int SearchOrder[]);//This function is used to generate individuals for initialisation

int randf(void)
{//just pick up a previously generated random integer
	if(random_i>999999)
	{
		random_i=0;
	}
	else
	{
		random_i++;
	}
	return random[random_i];
}
void makerand(void)
{
	int i,j;
	time_t t;
	srand((unsigned)time(&t));//set the random number seads
	for(i=0; i<1000000; i++)//make random numbers
	{
		random[i]=rand();
	}
	for (i = 0; i < FENCESIZE; i++)//make random fences
	{
		for(j = 0; j < 32; j++)
		{
			if(randf()%10 > 4)ADD_ELEMENT(& Fence[i],j);//I wish to get 50% 1s and 50% 0s in a 32 bit fence
			else DEL_ELEMENT(& Fence[i],j);
		}
	}
}
int Set_Size(int Set[])//to count the number of ones in the set.
{
	int i,d = 0;
	for (i = 0; i <  nCprEdg ;i++)
		d += POP_COUNT(Set[i]);
	return d;
}
int Read_Graph(int G[NMAX][MMAX])// to read the graph into the memory
{
	int i,j,u,d;
	nVerNum = nCprEdg = 0;
	scanf("%d",&nVerNum);
	if (nVerNum == 0)exit(0);
	nCprEdg = (nVerNum + 31) / 32;
	for(i = 0; i <  nVerNum; i++)
	{
		for(j = 0; j < nCprEdg; j++)
		{
			G[i][j] = 0;
		}
	}
	for(i = 0; i < nVerNum; i++)
	{
		scanf("%d",&d);
		for(j = 0; j < d; j++)
		{
			scanf("%d",&u);
			ADD_ELEMENT(G[i], u);
		}
		ADD_ELEMENT(G[i], i);
	}
	return 1;
}
int Print_Set(int n, int Set[], int j)//print the current set
{
	int i;
	for (i = 0; i < n; i++)
		if(IS_ELEMENT(Set,i) && i != j)printf(" %d", i);
	printf("\n");
	return 0;
}
int Can_Be_Blue(int G[NMAX][MMAX], int Blue[])//if the new added point can be blue(fixed it out of dominating set)return 1
{
	int i,j,d;
	for (i = 0; i < nVerNum; i++)
	// if any of the vertices can not be dominated by the time, the new added vertex must be removed from the blue set 
	{
		for (d = j = 0; j < nCprEdg; j++)
		// if all the neighbors(including itself) are 0(no connection) except those are in the blue set d will be 0.
			if(G[i][j] & ~Blue[j])
			// if there is still some neighbor which is not in the blue set, or if itself is not in the blue set, d will be 1
			{
				d = 1;
				break;
			}
		if(!d)return 0;
		// if d is 0, that means some vertex can not be dominated already, return 0.
	}
	return 1;// else return 1.
}
int Cross_Over(int Generation[CROSSOVERSIZE][NMAX])
{
	int nPos, i, j, k;
	nPos = GROUPSIZE;
	for (i = 1; i < GROUPSIZE + 1; i++)//the one has the highest fitness value will be able to crossover with all the other individuals
	{								//the one has the worst fitness value will be able to crossover with only the highest one.
		for(j = 0; j < i; j++)		//and the # of individuals that one can be mate with is determined by its ranking of the fitness value.
		{
			for (k = 0; k < nVerNum; k++)
			{
				if(IS_ELEMENT(&Fence[fence_i],(k % 32)))Generation[nPos][k] = Generation[i][k];
				else Generation[nPos][k] = Generation[j][k];//in a fence, 1 means this part is from the father, and 0 means this part is from the mother, 
			}
			nPos++;
		}
	}
	if(fence_i >= FENCESIZE)fence_i = 0;//pick up another fence.
	else fence_i++;
	return 0;
}
int Mutation(int G[NMAX][MMAX], int Generation[CROSSOVERSIZE][NMAX])
{
	int i,j,nMutPos,nMutGene,nCount;
	for(i = GROUPSIZE; i < CROSSOVERSIZE; i++)//the newly generated individuals will have some chance to mutate on one bit.
	{
		nMutPos = randf() % nVerNum;//randomly select the mutate position
		nMutGene = randf() % Set_Size(G[nMutPos]);//randomly select the mutate gene, so that It can still be dominated after mutation
		for(nCount = -1, j = 0; j < nVerNum; j++)//find that gene
		{
			if(IS_ELEMENT(G[nMutPos],j))nCount++;
			if(nCount == nMutGene)
			{
				Generation[i][nMutPos] = j;//mutate !
				break;
			}
		}
	}
	return 0;
}
int BFS_Array(int BfsArray[NMAX], int G[NMAX][MMAX], int Red[MMAX], int Root)
{//Using BFS search to generate a SearchOrder, and our search behavior will follow this order.
	int i,j,nParent,nPointer,FlagExist;
	//nParent:integer index for indicating the current parent position in the BFQ.
	//nPointer:integer index on the end of the BFQ
	//FlagExist:Flag to show if the vertex was been visited.
	nParent = nPointer = 0;
	BfsArray[nPointer++] = Root;//This is BFQ
	while(nParent < nPointer + 1)
	{
		for (i = 0; i < nVerNum; i++)
		{
			if(IS_ELEMENT(G[nParent], i) && !IS_ELEMENT(Red, i))//I traverse the parents' neighbors, and add them in to the BFQ.
			{//However, if some of the neighbor is already in the Red set, that means the neighbor has been fixed in dominating set for some reason, 
			//so I should ignore it when I are searching the min dom set.If it is not in our BFQ(will be the SearchOrder), the algorithm
			// will leave it alone.
				for(FlagExist = j = 0; j < nPointer; j++) if(BfsArray[j] == i)//before adding the neighbor into the BFQ, check if it is already in.
				{
					FlagExist = 1;
					break;
				}
				if(!FlagExist)BfsArray[nPointer++] = i;//if it is not in the current BFQ, Add it in!
			}
		}
		nParent++;
	}
	return 0;
}
int BinaryG_to_Generaiton(int G[NMAX][MMAX], int BinaryGen[], int Generation[])
{//Encode the domset from Compressed Adj Matrix to Genetic Chromosome Data Structure
	int i,j,Record;
	int Temp[MMAX];
	//Actually, we code each vertex with the vertex # of its parent ( Here I define the parent as the vertex which is dominating the child vertex, 
	//if one vertex is being dominated by mutiple vertices in the current dominating set, 
	//I will choose the max degree vertex as its dominating vertex to break tie.)
	for(i = 0; i < nVerNum; i++)
	{
		for(j = 0; j < nCprEdg; j++)Temp[j] = G[i][j] & BinaryGen[j];
		for(Record = -1,j = 0; j < nVerNum; j++)
		{
			if(IS_ELEMENT(Temp, j))
			{
				if(Record == -1)Record = j;
				else
					if(Set_Size(G[j]) > Set_Size(G[Record]))Record = j;//choose the maximum parent
			}
		}
		Generation[i] = Record;
	}
	return 0;
}
int Survival_Of_The_Fittest(int G[NMAX][MMAX], int Generation[CROSSOVERSIZE][NMAX])
{
	int i,j,k,SwampTemp,BubbleFlag;
	int BinaryGen[MMAX],RankTable[CROSSOVERSIZE];
//Begin Translate the Chromosome into dominating set, and evaluate the fitness value by calculatiin the set size
	for (j = 0; j < MMAX; j++)BinaryGen[j] = 0;
	for (i = 0; i < CROSSOVERSIZE; i++)
	{
		for(j = 0; j < nVerNum; j++)ADD_ELEMENT(BinaryGen, Generation[i][j]);
		RankTable[i] =  Set_Size(BinaryGen);
		for(j = 0; j < nCprEdg; j++)BinaryGen[j] = 0;
	}
//End Translation and evaluation-------------------------------------------------

//Begin----------Bubble Sort-----------------------------------------------------
	for(i = 0; i < CROSSOVERSIZE; i++)
	{//Bubble sort to rank the fitness value, as well as the position of the individuals.
		for(j = BubbleFlag = 0; j < CROSSOVERSIZE - 1; j++)
		{
			if(RankTable[j] > RankTable[j + 1])
			{
				SwampTemp = RankTable[j];
				RankTable[j] = RankTable[j + 1];
				RankTable[j + 1] = SwampTemp;
				for(k = 0; k < nVerNum; k++)
				{
					SwampTemp = Generation[j][k];
					Generation[j][k] = Generation[j + 1][k];
					Generation[j + 1][k] = SwampTemp;
				}
				BubbleFlag = 1;
			}
		}
		if(!BubbleFlag)break;
	}
//End----------Bubble Sort--------------------------------------------------------
	return RankTable[0];//return the best fitness value so far.
}
int Initialise(int G[NMAX][MMAX], int Generation[CROSSOVERSIZE][NMAX], int DomtedSet[])
{
	int i, j, k, BubbleFlag, SwampTemp, Order[2][NMAX];
	int BinaryGen[MMAX];
	int Red[MMAX];
	int Blue[MMAX];
	int SearchOrder[NMAX];
	for(i = 0; i < nVerNum; i++)
	{
		Order[0][i] = i;
		Order[1][i] = Set_Size(G[i]);
	}
	for(i = 0; i < nVerNum; i++)
	//rank the vertices according to the degree.
	{
		for(j = BubbleFlag = 0; j < nVerNum - 1; j++)
		{
			if(Order[1][j] < Order[1][j + 1])
			//bubble rank
			{
				SwampTemp = Order[1][j];
				Order[1][j] = Order[1][j + 1];
				Order[1][j + 1] = SwampTemp;
				SwampTemp = Order[0][j];
				Order[0][j] = Order[0][j + 1];
				Order[0][j + 1] = SwampTemp;
				BubbleFlag = 1;
			}
		}
		if(!BubbleFlag)break;
	}
	for (i = 0; i < CROSSOVERSIZE; i++)
	{
		for (j = 0; j < MMAX; j++)//Reset
		{
			Red[j] = 0;
			Blue[j] = 0;
			DomtedSet[j] = 0;
			BinaryGen[j] = 0;
		}
		for(j = 0; j < nVerNum; j++)
		{
			SearchOrder[j] = 0;
			if(Order[1][j] == 2)
			//if the vertex just has degree 1(2 ones in compressed adjacency matrix), then its neighbor must be part of the dominating set.
			{
				for(k = 0; k < nVerNum; k++)if(IS_ELEMENT(G[j],k) && j != k) ADD_ELEMENT(Red,k);
			}
			else
			{
				if(Order[1][j] == 1) ADD_ELEMENT(Red,j);// if the vertex is separated, then it must be part of the dominating set.
			}
		}
		BFS_Array(SearchOrder, G, Red, Order[0][(i % nVerNum)]);
		Generate_Dom_Set(G, Red, Blue, BinaryGen, DomtedSet, 0, SearchOrder);
		BinaryG_to_Generaiton(G, BinaryGen, Generation[i]);
	}
	return 0;
}
void Generate_Dom_Set(int G[NMAX][MMAX], int Red[], int Blue[], int DomtingSet[], int DomtedSet[], int nDepthIndex,int SearchOrder[NMAX])
{
	// recursive function for finding the minimum dominating set
	int i, nUnDom, nDom;
	ADD_ELEMENT(Blue,SearchOrder[nDepthIndex]);//firstly, make the current vertex blue(out of dominating set).
	if(Can_Be_Blue(G,Blue))
	{
		Generate_Dom_Set(G, Red, Blue, DomtingSet, DomtedSet, nDepthIndex + 1, SearchOrder);
		return;
	}
	DEL_ELEMENT(Blue, SearchOrder[nDepthIndex]);
	ADD_ELEMENT(Red, SearchOrder[nDepthIndex]);
	//After trying blue, I make it Red(part of dominating set).
	for(i = 0; i < nCprEdg; i++) DomtedSet[i] = DomtedSet[i] | G[SearchOrder[nDepthIndex]][i];
	//After making the vertex a part of dominating set, I want to update the dominated set by now.
	nUnDom =nVerNum - Set_Size(DomtedSet);
	//check if all the vertices are dominated.
	nDom = Set_Size(Red);//calculate the size of the current dominating set.
	if(nUnDom == 0)// if the number of undominated vertices is zero, I will get a dominating set.
	{
		for(i = 0; i < nCprEdg; i++)DomtingSet[i] = Red[i];//copy set.
		return;
	}
	//if the program runs here, it means the number of undominated vertices is NOT zero.
	Generate_Dom_Set(G, Red, Blue, DomtingSet, DomtedSet, nDepthIndex + 1, SearchOrder);
	return;
}
int main(void)
{
	int nCount, i, Result;
	int G[NMAX][MMAX];
	int Generation[CROSSOVERSIZE][NMAX];
	int DomtedSet[MMAX];
	int MinDomtingSet[MMAX];
	nCount = 1;
	makerand();
	while(Read_Graph(G))
	{
		printf("%d %d\n", nCount,nVerNum);//print
		Initialise(G,Generation,DomtedSet);//initialization
		for (i = 0; i < MAX_LOOP * nVerNum; i++)
		{
			Result = Survival_Of_The_Fittest(G, Generation);
			Cross_Over(Generation);
			Mutation(G,Generation);
		}
		printf("2 %d  ",Result);
		for(i = 0; i < MMAX; i++)MinDomtingSet[i] = 0;
		for(i = 0; i < nVerNum; i++)ADD_ELEMENT(MinDomtingSet, Generation[0][i]);
		Print_Set(nVerNum,MinDomtingSet,-1);
		nCount++;
	}
	return 0;
}
//-------------------------------------------------by Ye Liang at 18:38 on 16th Dec. 2014-------------------------------------------------------------

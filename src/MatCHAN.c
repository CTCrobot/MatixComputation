//======================================================//
//                       _oo0oo_                        //
//                      o8888888o                       //
//                      88" . "88                       //
//                      (| -_- |)                       //
//                      0\  =  /0                       //
//                    ___/`---'\___                     //
//                  .' \\|     |// '.                   //
//                 / \\|||  :  |||// \                  //
//                / _||||| -:- |||||- \                 //
//               |   | \\\  -  /// |   |                //
//               | \_|  ''\---/''  |_/ |                //
//               \  .-\__  '-'  ___/-. /                //
//             ___'. .'  /--.--\  `. .'___              //
//          ."" '<  `.___\_<|>_/___.' >' "".            //
//         | | :  `- \`.;`\ _ /`;.`/ - ` : | |          //
//         \  \ `_.   \_ __\ /__ _/   .-` /  /          //
//     =====`-.____`.___ \_____/___.-`___.-'=====       //
//                       `=---='                        //
//                                                      //
//                                                      //
//     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~      //
//                                                      //
//               佛祖保佑         永无BUG                //
//======================================================//

#include "MatCHAN.h"

 //==================================================================
//函数名：    mc_fnvMatrixPrint
//功能：      矩阵打印函数
//==================================================================
void mc_fnvMatrixPrint(double **dInMatrix, int nRow, int nColumn) {
	//printf("\n[");
	for (int i = 0; i < nRow; i++) {
		for (int k = 0; k < nColumn; k++)
		{
			printf(" %4.4f ", *((double *)dInMatrix + i * nColumn + k));
		}
		printf("\n");
	}
	//printf("\b\b\b\b]\n");
}

 //==================================================================
//函数名：    mc_fnvMatrixPrint
//功能：      向量打印函数
//==================================================================
void mc_fnvVectorPrint(double *dInVector, int nNumber) {
	for (int i = 0; i < nNumber; i++) 
		printf(" %4.4f ", *(dInVector + i));
    printf("\n");
}

//===================================================================
//===================================================================
//对于矩阵相关运算 需要提前设定一个可能的最大矩阵帮助求解 脑子现在转不过来 malloc用的不熟先这么着吧
//使用预设的数组 很方便
double dSubMatrix[MAX_DIM][MAX_DIM] = { 0.0 };
 //==================================================================
 //函数名：    mc_fnfMatrixDeterminant
 //功能：      求矩阵的行列式 采用递归的方法
 //==================================================================
double mc_fnfMatrixDeterminant(double dInMatrix[MAX_DIM][MAX_DIM], int nRow, int nColumn) {
	/* 	function that calculate the determinant of a matrix
		using the Laplace algorithm
		double **m pointer to the matrix
		int n  dimension of a matrix
	*/
	double dFunMatrix[MAX_DIM][MAX_DIM];
	double determinante = 0.0;
	int row, column, i, j;

	// 将矩阵复制一次
	for (int i = 0; i < nRow; i++)
		for (int k = 0; k < nColumn; k++)
			dFunMatrix[i][k] = dInMatrix[i][k];

	if (nRow >= nColumn) {
		/*determiant of a 1x1 matrix*/
		if (nRow == 1 && nColumn == 1) {
			determinante = dFunMatrix[0][0];
		}
		else if (nRow == 2 && nColumn == 2) {
			/*determinant of a 2x2 matrix*/
			determinante = dFunMatrix[0][0] * dFunMatrix[1][1] - dFunMatrix[0][1] * dFunMatrix[1][0];
		}
		else {
			/*determinant of a matrix upper 2x2 dimension*/
			for (row = 0; row < nRow; row++) {
				/*	create a submatrix that will be used in input
					in fnfMatrixDeterminant() function
				*/
				for (i = 0; i < nRow - 1; i++) {
					for (j = 0; j < nColumn - 1; j++) {
						int sub_row = (i < row ? i : i + 1);
						/*	apply the algorithm for the first column*/
						int sub_col = j + 1;
						dSubMatrix[i][j] = dInMatrix[sub_row][sub_col];
					}
				}
				if (row % 2 == 0) {
					determinante += dFunMatrix[row][0] * mc_fnfMatrixDeterminant(dSubMatrix, nRow - 1, nColumn - 1);
				}
				else {
					determinante -= dFunMatrix[row][0] * mc_fnfMatrixDeterminant(dSubMatrix, nRow - 1, nColumn - 1);
				}
			}
		}
	}
	else {
		/*determinant of a matrix upper 2x2 dimension*/
		for (column = 0; column < nColumn; column++) {
			for (i = 0; i < nRow - 1; i++) {
				for (j = 0; j < nColumn - 1; j++) {
					int sub_row = i + 1;
					int sub_col = (j < column ? j : j + 1);
					dSubMatrix[i][j] = dInMatrix[sub_row][sub_col];
				}
			}
			if (column % 2 == 0) {
				determinante += dFunMatrix[0][column] * mc_fnfMatrixDeterminant(dSubMatrix, nRow - 1, nColumn - 1);
			}
			else {
				determinante -= dFunMatrix[0][column] * mc_fnfMatrixDeterminant(dSubMatrix, nRow - 1, nColumn - 1);
			}
		}
	}
	
	return determinante;
}

//==================================================================
//函数名：    mc_fnvMatrixTransport
//功能：      矩阵转置 可以转置任意大小的矩阵 dDet 变量是矩阵的行列式 只求转置时赋值1即可
//说明：      nSizeMax 这个参数非常重要  如果是提前设定了矩阵大小 但是放的实际容量有小于这个量就特别注意了
//
//==================================================================
void mc_fnvMatrixTransport(double **dInRotMatrix, double **dOutRotMatrix,int nSizeMax, int nRow, int nColuw, double dDet) {
	for (int i = 0; i < nColuw; i++) {
		for (int k = 0; k < nRow; k++) {
			*((double*)dOutRotMatrix+i*nRow + k) = (*((double*)dInRotMatrix + k * nSizeMax + i))/ dDet;
		}
	}
}

//==================================================================
//函数名：    mc_fnvMatrixInverse
//功能：      方阵矩阵求逆 由于知道自己传入的矩阵是好的 所以没有检测错误
//==================================================================
void mc_fnvMatrixInverse(double **dInRotMatrix, double **dOutRotMatrix, int n) {
	/*	function thet calculate the cofactors of a matrix
		double **m pointer to the matrix
		double d[MAX_DIM][MAX_DIM] temporary matrix
		int n dimension
		double det determinant
	*/
	double c[MAX_DIM][MAX_DIM] = {0.0}, dTemp[MAX_DIM][MAX_DIM];
	double b[MAX_DIM][MAX_DIM];
	int l, h, r, k, i, j;

	// 将矩阵复制一次
	for (int i = 0; i < n; i++)
		for (int k = 0; k < n; k++)
			dTemp[i][k] = (*((double*)dInRotMatrix + i * 4 + k));

	double dDet = mc_fnfMatrixDeterminant(dTemp, n,n);

	for (h = 0; h < n; h++) {
		for (l = 0; l < n; l++) {
			r = 0;
			k = 0;
			for (i = 0; i < n; i++)
				for (j = 0; j < n; j++)
					if (i != h && j != l) {
						b[r][k] = (*((double*)dInRotMatrix + i*n + j));
						if (k < (n - 2))
							k++;
						else {
							k = 0;
							r++;
						}
					}

			c[h][l] = pow(-1, (h + l))*mc_fnfMatrixDeterminant(b, (n - 1), (n - 1));
		}
	}
	mc_fnvMatrixTransport(c, dOutRotMatrix, MAX_DIM, n,n, dDet);
}
//===================== 求逆矩阵用的相关函数 =========================
//==================================================================

//==================================================================
//函数名：    fndVectorNorm
//功能：      向量求模
//==================================================================
double mc_fndVectorNorm(double *dVector,int nNumber){ 
	double dNorm = 0.0;
	for(int i = 0; i < nNumber; i++)
		dNorm = dNorm + pow(*(dVector+i),2);

	dNorm = sqrt(dNorm);
	return dNorm;
}

//==================================================================
//函数名：    mc_fndMatrixMultp
//功能：      向量矩阵的相乘 基于斯特拉森方法适用于维数大于8阶的情况 低维就传统算法吧
//说明：      之所以第二个矩阵不用传进来行数 是因为它必然等于矩阵1的列数
//==================================================================
void mc_fnvMatrixMultp(double **dMatrix1,double **dMatrix2, double **dOutMatrix, int nRow1, int nColumn1, int nColumn2){ 
	for(int i = 0; i<nRow1; i++){
        for(int k = 0; k<nColumn2; k++){
            for(int j = 0; j<nColumn1; j++){
                (*((double*)dOutMatrix + i * nColumn2 + k)) += (*((double*)dMatrix1 + i * nColumn1 + j))*(*((double*)dMatrix2 + j * nColumn2 + k));
            }
        }
    }
}

//==================================================================
//函数名：    mc_fnvMatrixVectorMultp
//功能：      矩阵与向量相乘
//说明：      
//==================================================================
void mc_fnvMatrixVectorMultp(double **dMatrix, double *dVector, double *dOutVector, int nRow, int nColumn) {
	for (int i = 0; i < nRow; i++) {
		for (int k = 0; k < nColumn; k++) {
			(*(dOutVector + i)) += (*((double*)dMatrix + i * nColumn + k))*(*(dVector + k));
		}
	}
}

//==================================================================
//函数名：    mc_fnvEularToRotMatrix
//功能：      将xyz的欧拉角转换为旋转矩阵 右手系
//==================================================================
void mc_fnvEularToRotMatrix(double dEl[3], double dOutRotMatrix[3][3]) { // 这个函数可以在后期设置成求解IK的约束 比如有些关节只能在一定范围内活动
	double qx = dEl[0], qy = dEl[1], qz = dEl[2];

	dOutRotMatrix[2][2] = cos(qy)*cos(qz);
	dOutRotMatrix[2][1] = cos(qy)*sin(qz);
	dOutRotMatrix[2][0] = -sin(qy);

	dOutRotMatrix[1][2] = -cos(qx)*sin(qz) + cos(qz)*sin(qx)*sin(qy);
	dOutRotMatrix[1][1] = cos(qx)*cos(qz) + sin(qz)*sin(qx)*sin(qy);
	dOutRotMatrix[1][0] = cos(qy)*sin(qx);

	dOutRotMatrix[0][2] = sin(qx)*sin(qz) + cos(qz)*cos(qx)*sin(qy);
	dOutRotMatrix[0][1] = -sin(qx)*cos(qz) + sin(qz)*cos(qx)*sin(qy);
	dOutRotMatrix[0][0] = cos(qx)*cos(qy);
}

//==================================================================
//函数名：    mc_fnvRotMatrixToOmega
//功能：      将旋转矩阵转化为ω
//说明：      
//==================================================================
void mc_fnvRotMatrixToOmega(double dMatrix[3][3], double dOutOmega[3]) {
	double dE1[3] = { dMatrix[2][1] - dMatrix[1][2],dMatrix[0][2] - dMatrix[2][0],dMatrix[1][0] - dMatrix[0][1] };
	double dNormE1 = mc_fndVectorNorm(dE1, 3);
	if (dNormE1 > EPS_MC)
		for (int i = 0; i < 3; i++)
			dOutOmega[i] = atan2(dNormE1, dMatrix[0][0] + dMatrix[1][1] + dMatrix[2][2] - 1) / dNormE1 * dE1[i];
	else if (dMatrix[0][0] > 0.0 && dMatrix[1][1] > 0.0 && dMatrix[2][2] > 0.0)
		for (int i = 0; i < 3; i++)
			dOutOmega[i] = 0.0;
	else
		for (int i = 0; i < 3; i++)
			dOutOmega[i] = (dMatrix[i][i] + 1)*3.1415926 / 2;
}

//==================================================================
//函数名：    mc_fnvDiagMatrix
//功能：      将所给向量转化为对角矩阵形式
//说明：      dOutMatrix需要提前声明为全零矩阵
//==================================================================
void mc_fnvDiagMatrix(double *dInVector, double **dOutMatrix, int nNumber) {
	for (int i = 0; i < nNumber; i++)
	{
		(*((double*)dOutMatrix + i * nNumber + i)) = *(dInVector + i);
	}
}

//==================================================================
//函数名：    mc_fnvIdentityMatrix
//功能：      产生一个单位矩阵*dMultNumber+dPlusNumber[i]
//说明：      dOutMatrix需要提前声明为全零矩阵
//==================================================================
void mc_fnvIdentityMatrix( double **dOutMatrix, int nNumber, double dMultNumber, double *dPlusNumber) {
	for (int i = 0; i < nNumber; i++)
	{   
        for(int k = 0; k<nNumber;k++){
            if(i==k)
		        (*((double*)dOutMatrix + i * nNumber + i)) = dMultNumber+(*(dPlusNumber+i));
            else
                (*((double*)dOutMatrix + i * nNumber + i)) = 0.0;
        }
	}
}

//==================================================================
//函数名：    mc_fndVectorCopy
//功能：      将向量拷贝返回同样的维度的向量
//说明：      
//==================================================================
double* mc_fndVectorCopy(double *dInVector, int nNumber) {
	double *dOutVector = (double *)malloc(sizeof(double)*nNumber);
	for (int i = 0; i < nNumber; i++)
	{
		*(dOutVector + i) = *(dInVector + i);
	}
	return dOutVector;
}

//==================================================================
//函数名：    mc_fndMatixCopy
//功能：      将矩阵拷贝返回同样的维度的矩阵
//说明：      
//==================================================================
double** mc_fndMatixCopy(double **dInMatrix, int nRow, int nColumn) {

	double **dOutMatrix = (double **)malloc(sizeof(double*)*nRow);
	for(int i = 0; i<nRow; i++)
		dOutMatrix[i] = (double *)malloc(sizeof(double)*nColumn);

	for (int i = 0; i < nRow; i++)
	{
		for (int k = 0; k < nColumn; k++)
		{
			*((double*)dOutMatrix+i*nColumn+k) = *((double*)dInMatrix + i * nColumn + k);
		}
	}
	return dOutMatrix;
}

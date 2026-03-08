#include<iostream>
#include<random>
#include<cmath>
#include<algorithm>
#include<set>
#include<unordered_set>

using namespace std;

typedef int F2;
typedef int F3;

const int verbose = true;

// the length of key and input
const int NN = 32; // 64
                   //
// security: conservative NN / 2.5; agressive: NN / 2
const int lambda = NN / 2;

// n/2
const int HalfNN = NN / 2; // 32

// output length
const int TT = int( lambda / log2(3) ); // 32 / 1.5 = 24

// for n = 2^m
const int RR = ( TT + 1 ) / 2; // 12 

// for n = 4m
//const int RR = NN / 4; // 12 

//const int RR = lambda / 2; // 12 
                       //
                       //
/* utility functions */
template< typename T >
void printMatrix( T ** M, int NRow, int NCol )
{
    for ( int i = 0; i < NRow; i++ )
    {
        for ( int j = 0; j < NCol; j++ )
            cout << static_cast<int> ( M[i][j] );
        cout << endl;
    }
}

inline F3 mod3(F3 x) 
{ 
    return (x % 3 + 3) % 3; 
}

inline F3 inverse_mod3(F3 x) 
{ 
    return x == 1 ? 1 : 2; 
} 

void generateCombinations( vector<vector<F3>>& res, vector<F3>& current, int depth, int k) 
{
    if (depth == k) 
    {
        res.push_back(current);
        return;
    }
    for (F3 val = 0; val < 3; ++val) 
    {
        current[depth] = val;
        generateCombinations(res, current, depth + 1, k);
    }
}

bool EQ( int * A, int * B, int N )
{
    for ( int i = 0; i < N; i++ )
        if ( A[i] != B[i] )
            return false;
    return true;
}

bool EQ( int * A, const vector<int> & B, int N )
{
    for ( int i = 0; i < N; i++ )
        if ( A[i]  != B[i] )
            return false;
    return true;
}

bool EQ( const vector<int> & A, int * B, int N )
{
    for ( int i = 0; i < N; i++ )
        if ( A[i]  != B[i] )
            return false;
    return true;
}

bool EQ( const vector<int> & A, const vector<int> & B, int N )
{
    for ( int i = 0; i < N; i++ )
        if ( A[i]  != B[i] )
            return false;
    return true;
}


/*
inline ostream & operator<< ( ostream & os, int t )  
{
    os << static_cast<int> ( t );
    return os;
}
*/

void printVec( int * t, int N )
{
    for ( int i = 0; i < N; i++ )
        cout << t[i];
    cout << endl;
}

void printVec( vector<int> & t, int N )
{
    for ( int i = 0; i < N; i++ )
        cout << t[i];
    cout << endl;
}


/* the wPRF implementation */

void MatrixMulVector2( F2 ** K, F2 * x, F2 * y, int N = NN )
{
    for ( int i = 0; i < N; i++ )
    {
        y[i] = 0;
        for ( int j = 0; j < N; j++ )
            y[i] = ( y[i] + K[i][j] * x[j] ) % 2;
    }
}

void MatrixMulVector3( F3 ** B, F3 * x, F3 * y, int T = TT, int N = NN )
{
    for ( int i = 0; i < T; i++ )
    {
        y[i] = 0;

        for ( int j = 0; j < N; j++ )
            y[i] = ( y[i] + B[i][j] * x[j] ) % 3;
    }
}

void wPRF( F2 ** K, F3 ** B, F2 * x, F3 * y )
{
    F2 * xx = new F2[NN];

    MatrixMulVector2( K, x, xx );
    MatrixMulVector3( B, xx, y );

    delete [] xx;
}

/* the wPRF implementation finished */

/* offline phase */

int rankF3( F3** B, int TT, int NN) 
{
    F3** mat = new F3*[TT];
    for (int i = 0; i < TT; ++i ) 
    {
        mat[i] = new F3[NN];
        for (int j = 0; j < NN; ++j) 
            mat[i][j] = B[i][j]; 
    }

    int rank = 0;
    for (int col = 0; col < NN && rank < TT; ++col) 
    {
        int pivot_row = -1;
        for (int i = rank; i < TT; ++i) 
        {
            if (mat[i][col] != 0) 
            {
                pivot_row = i;
                break;
            }
        }
        if (pivot_row == -1) 
            continue; 

        //int* temp = mat[rank];
        //mat[rank] = mat[pivot_row];
        //mat[pivot_row] = temp;

        swap( mat[rank], mat[ pivot_row ] );

        int inv = inverse_mod3( mat[rank][col] );
        for (int j = 0; j < NN; ++j) 
            mat[rank][j] = mod3(mat[rank][j] * inv);

        for (int i = rank + 1; i < TT; ++i) 
        {
            if ( mat[i][col] == 0 ) 
                continue;
            int factor = mod3( -mat[i][col] );
            for (int j = col; j < NN; ++j) 
                mat[i][j] = mod3( mat[i][j] + ( factor * mat[rank][j] ) );
        }
        ++rank;
    }

    for (int i = 0; i < TT; ++i) {
        delete[] mat[i];
    }
    delete[] mat;

    return rank;
}

int rankF2(F2** KEY, int NN) 
{
    F2** mat = new F2*[NN];
    for (int i = 0; i < NN; ++i) 
    {
        mat[i] = new F2[NN];
        for (int j = 0; j < NN; ++j) 
            mat[i][j] = KEY[i][j];
    }

    int rank = 0;
    for (int col = 0; col < NN && rank < NN; ++col) 
    {
        int pivot_row = -1;
        for (int i = rank; i < NN; ++i) 
        {
            if (mat[i][col] != 0) 
            {
                pivot_row = i;
                break;
            }
        }
        if (pivot_row == -1) 
            continue;

        swap( mat[rank], mat[ pivot_row ] );

        for (int i = 0; i < NN; ++i) 
        {
            if (i == rank || mat[i][col] == 0) 
                continue;
            for (int j = col; j < NN; ++j) 
                mat[i][j] = ( mat[i][j] + mat[rank][j] ) % 2;
        }
        ++rank;
    }

    for (int i = 0; i < NN; ++i) 
        delete[] mat[i];
    delete[] mat;

    return rank;
}

int Gauss_F2( F2 ** mat, int NRow, int NCol )   
{
    for ( int col = 0; col < NCol; col++ )
    {
        int currentRow = col;

        // current row, find the first pivot
        for ( int row = currentRow; row < NRow; row ++ )
        {
            if ( mat[row][col] == 0 )
                continue;
            else
                //swapRow( mat, NCol, currentRow, row );
                swap( mat[ currentRow ], mat[ row ] );
        }

        for ( int row = currentRow + 1; row < NRow; row ++ )
        {
            if ( mat[row][col] == 1 )
                for ( int j = 0; j < NCol; j++ )
                    mat[row][j] = ( mat[row][j] + mat[currentRow][j] ) % 2;
        }
    }

    // eliminate the above elements of the pivot
    for ( int col = 0; col < NCol; col++ )
    {
        if ( mat[col][col] == 1 )
        {
            for ( int i = 0; i < col; i++ ) // for all rows above pivot
            {
                if ( mat[i][col] == 1 )
                    for ( int j = 0; j < NCol; j++ )
                        mat[i][j] = ( mat[i][j] + mat[col][j] ) % 2; 
            }
        }
    }

    int rank = 0;
    for ( int i = 0; i < NCol; i++ )
        if ( mat[i][i] == 1 )
            rank += 1;

    return rank;
}

//void genGuesseDiffs ( F2 ** mat, int NRow, int NCol, F2 ** diff ) 
vector<vector<F2> > genGuesseDiffs ( F2 ** mat, int NRow, int NCol ) 
{
    //printMatrix( mat, NRow, NCol );
    //getchar();

    vector< vector<F2> > diff;

    F2 * X = new F2 [ NCol ]; 

    for ( int x = 0; x < ( 1 << NCol ); x++ )
    {
        for ( int i = 0; i < NCol; i++ )
            X[i] = x >> ( NCol - 1 - i ) & 0x1; 

        vector<F2> vec;

        for ( int i = 0; i < NRow; i++ )
        {
            F2 temp = 0;
            for ( int j = 0; j < NCol; j++ )
                if ( mat[i][j] == 1 )
                    temp = ( temp + X[j] ) % 2;   
            vec.push_back( temp );
        }

       // for ( int i = 0; i < HalfNN; i++ )
       //     cout << vec[i];
       // cout << endl;
       // getchar();

        diff.push_back( vec );
    }

    delete [] X;

    return diff;
}
/* offline phase finished */

/* solve a F3 linear equation systems  */

// T rows and N cols
vector<vector<F3>> solveF3(F3** A, F3* b, int T, int N) 
{
    vector<vector<F3>> aug(T, vector<F3>(N + 1));

    // init the augmented matrix
    for (int i = 0; i < T; ++i) 
    {
        for (int j = 0; j < N; ++j) 
            aug[i][j] = mod3(A[i][j]);
        aug[i][N] = mod3(b[i]);
    }

    vector<int> pivot_cols;
    int rank = 0;

    for (int col = 0; col < N; ++col) 
    {
        int pivot = rank;
        while ( ( pivot < T ) && ( aug[pivot][col] == 0 ) ) 
            ++pivot;
        if (pivot == T) 
            continue; // free variable
                      //

        swap( aug[rank], aug[pivot] );


        int inv = inverse_mod3(aug[rank][col]);


        for (int j = 0; j <= N; ++j)
            aug[rank][j] = mod3( aug[rank][j] * inv );


        for (int i = 0; i < T; ++i) 
        {
            if ( ( i != rank ) && ( aug[i][col] != 0 ) ) 
            {
                int factor = mod3(-aug[i][col]);
                for (int j = col; j <= N; ++j)
                    aug[i][j] = mod3(aug[i][j] + factor * aug[rank][j]);
            }
        }
        pivot_cols.push_back(col);
        ++rank;
    }

    for (int i = rank; i < T; ++i)
        if (aug[i][T] != 0) 
            return {}; // no solution

    unordered_set<int> pivot_set(pivot_cols.begin(), pivot_cols.end());

    vector<int> free_cols;
    for (int col = 0; col < N; ++col)
        if (!pivot_set.count(col)) 
            free_cols.push_back(col);

    vector<vector<F3>> solutions;

    int k = free_cols.size();
    vector<vector<F3>> free_combos;
    if ( k > 0 ) 
    {
        vector<F3> tmp(k, 0);
        generateCombinations(free_combos, tmp, 0, k);

        for (auto& combo : free_combos) 
        {
            vector<F3> sol(N, 0);
            for (int i = 0; i < k; ++i) 
                sol[free_cols[i]] = combo[i];

            for (int i = rank-1; i >= 0; --i) 
            {
                int pc = pivot_cols[i];
                F3 sum = 0;
                for (int j = pc+1; j < N; ++j)
                    sum = mod3(sum + mod3(aug[i][j] * sol[j]));
                sol[pc] = mod3(aug[i][N] - sum);
            }
            solutions.push_back(sol);
        }
    }
    else
    {
        vector<F3> sol(N, 0);
        for (int i = rank - 1; i >= 0; --i) {
            int pc = pivot_cols[i];
            F3 sum = 0;
            for (int j = pc + 1; j < N; ++j)
                sum = mod3(sum + mod3(aug[i][j] * sol[j]));
            sol[pc] = mod3(aug[i][N] - sum);
        }
        solutions.push_back(sol);
    }

    return solutions;
}


void generateCombinationsF2( vector<vector<F2>>& res, vector<F2>& current, int depth, int k) 
{
    if (depth == k) 
    {
        res.push_back(current);
        return;
    }
    for (F2 val = 0; val < 2; ++val) 
    {
        current[depth] = val;
        generateCombinationsF2(res, current, depth + 1, k);
    }
}


// M is the matrix, b is the b, T is the number of rows, N is the number of columns
vector<vector<F2>> solveF2(F2** M, F2* b, int T, int N) 
{
    vector<vector<F2>> aug(T, vector<F2>(N + 1));

    for (int i = 0; i < T; ++i) {
        for (int j = 0; j < N; ++j)
            aug[i][j] = M[i][j];

        aug[i][N] = b[i];
    }

    vector<int> pivot_cols; 
                        
    int rank = 0;

    for (int col = 0; col < N; ++col) 
    {
        int pivot = rank;

        while ( ( pivot < T ) && ( aug[pivot][col] == 0 ) )
            ++pivot;
        if (pivot == T)
            continue; 

        swap(aug[rank], aug[pivot]);

        for (int i = 0; i < T; ++i) 
            if ( ( i != rank ) && ( aug[i][col] != 0 ) ) 
                for (int j = col; j <= N; ++j)
                    aug[i][j] = aug[i][j] ^ aug[rank][j];

        pivot_cols.push_back(col);
        ++rank;
    }

    for (int i = rank; i < T; ++i) 
        if (aug[i][N] != 0)
            return {}; 

    unordered_set<int> pivot_set(pivot_cols.begin(), pivot_cols.end());
    vector<int> free_cols;
    for (int col = 0; col < N; ++col) {
        if (!pivot_set.count(col))
            free_cols.push_back(col);
    }

    vector<vector<F2>> solutions;
    int k = free_cols.size(); 

    if (k > 0) 
    {
        vector<vector<F2>> free_combos;
        vector<F2> tmp(k, 0);
        generateCombinationsF2(free_combos, tmp, 0, k);

        for (auto& combo : free_combos) {
            vector<F2> sol(N, 0);
            
            for (int i = 0; i < k; ++i)
                sol[free_cols[i]] = combo[i];
            for (int i = rank - 1; i >= 0; --i) {
                int pc = pivot_cols[i]; 
                F2 sum = 0;
                for (int j = pc + 1; j < N; ++j)
                    sum = sum ^ aug[i][j] * sol[j];
                sol[pc] = aug[i][N] ^ sum;
            }
            solutions.push_back(sol);
        }
    } else {
        vector<F2> sol(N, 0);
        for (int i = rank - 1; i >= 0; --i) {
            int pc = pivot_cols[i];
            F2 sum = 0;
            for (int j = pc + 1; j < N; ++j)
                sum = sum ^ aug[i][j] * sol[j];
            sol[pc] = aug[i][N] ^ sum;
        }
        solutions.push_back(sol);
    }

    return solutions;
}

int main()
{
    int S1 = 0;
    int S2 = 0;
    int W1_before = 0;
    int W1 = 0;
    int W2 = 0;
    int E = 0;
    int SUC = 0;

    random_device rd;
    mt19937 gen ( rd() );
    uniform_int_distribution<F2> dis(0, 1);
    uniform_int_distribution<F3> disF3(0, 2);

    float TESTNUM = 100;

    int COUNT = 0;

    cout << "The parameter NN = " << NN << ";" << "lambda = " << lambda << endl;

    cout << "Execute 100 times of tests...The data in the paper is the average value of the 100 tests..." << endl;

    while ( true )
    {
        if ( verbose )
        {
            cout << endl;
            cout << "Test time: " << COUNT << endl;
        }


        if ( verbose )
        {
            cout << "-----------------------------------------------------------" << endl;
            cout << "Init the wPRF " << endl;
            cout << "PRF parameters: N = " << NN << " T = " << TT  << " lambda = " << lambda << endl;
        }
    

        /* initialize the wPRF */ 
        /* generate the key */
        F2 ** KEY = new F2 * [NN];
        for ( int i = 0; i < NN; i++ ) 
            KEY[i] = new F2[NN];

        for ( int i = 0; i < NN; i++ )
            KEY[0][i]= dis( gen ); // init the key
                                   //
        for ( int i = 1; i < NN; i++ )
            for ( int j = 0; j < NN; j++ ) 
                KEY[i][ j ] = KEY[0][ ( i + j ) % NN ];

        //printMatrix( KEY, NN, NN );

        // compute the rank of KEY
        int rank_key = rankF2( KEY, NN );   

        if ( verbose )
        {
            if ( rank_key != NN )
            {
                cout << "The key rank is not full, run again" << endl;
                //return -1;
                continue;
            }
        }

        COUNT += 1;

        if ( verbose )
        {
            cout << "Init the key... Done! " << endl;
            for ( int i = 0; i < NN; i++ )
                cout << static_cast<int>( KEY[0][i] );
            cout << endl;
        }

        F3 ** B = new F3 * [TT];
        for ( int i = 0; i < TT; i++ )
            B[i] = new F3[ NN ];

        for ( int i = 0; i < TT; i++ )
            for ( int j = 0; j < NN; j++ )
                B[i][j] = disF3( gen );

        int rankB = rankF3( B, TT, NN );

        if ( verbose)
        {
            if ( rankB < TT )
            {
                cout << "B does not reach TT rank " << endl;
                return -2;
            }
        }

        if ( verbose)
        {
            cout << "Init the B matrix... Done! " << endl;
        }

        //printMatrix<F3> ( B, TT, NN );

        if ( verbose )
        {
            // attack parameter 
            cout << "-----------------------------------------------------------" << endl;
            cout << "Offline Phase " << endl; 

            cout << "Attack Rank R = " << RR << endl;

        }

        F2 ** mat = new F2 * [ HalfNN ];
        for ( int i = 0; i < HalfNN; i++ )
            mat[i] = new F2[ HalfNN ];

        F2 ** matT = new F2 * [HalfNN];
        for ( int i = 0; i < HalfNN; i++ )
            matT[i] = new F2[ RR ];

        F2 ** V = new F2 * [HalfNN];
        for ( int i = 0; i < HalfNN; i++ )
            V[i] = new F2[ HalfNN ];

        if ( verbose )
        {
            cout << "Search for an input with difference rank = " << RR << endl;
            cout << "The expected complexity is (Proposition 1) 2^{R-1-N/2} = 2^" << ( RR - 1 - HalfNN ) << endl;
        }

        F2 * Input1 = new F2[ NN ];
        F3 * Output1 = new F3[ TT ];

        F2 * Input2 = new F2[ NN ];
        F3 * Output2 = new F3[ TT ];

        F2 ** U1 = new F2 * [ HalfNN + RR ];
        for ( int i = 0; i < HalfNN + RR; i++ )
            U1[i] = new F2[ NN ];

        for ( int i = 0; i < HalfNN + RR; i++ )
            for ( int j = 0; j < NN; j++ )
                U1[i][j] = 0;

        for ( int i = 0; i < HalfNN + RR; i++ )
            U1[i][i] = 1;

        if ( verbose )
        {
            cout << "U1 is " << endl;
            printMatrix( U1, RR + HalfNN, NN );
        }
        

        //getchar();

        F2 ** U2 = new F2 * [ HalfNN + RR ];
        for ( int i = 0; i < HalfNN + RR; i++ )
            U2[i] = new F2[ NN ];

        long long query = 0;

        vector< vector<F2> > diff1; // diff space
        vector< vector<F2> > diff2; // diff space

        for ( query = 0; query < ( 1LL << 32 ); query++ )
        {
            // the input we are waiting for
            for ( int i = 0; i < NN; i++ )
                Input1[i] = dis( gen );

            for ( int i = 0; i < HalfNN; i++ )
                V[0][i] = ( Input1[i] + Input1[ i + HalfNN ] ) % 2;

            // the difference rotations
            for ( int row = 1; row < HalfNN; row++ )
                for ( int j = row; j < HalfNN + row; j++ )
                    V[row][j % HalfNN] = V[0][ j - row ];

            for ( int i = 0; i < HalfNN; i++ )
                for ( int j = 0; j < HalfNN; j++ )
                    mat[j][i] = V[i][j];

            int rank = Gauss_F2( mat, HalfNN, HalfNN );   

            if ( rank == RR )
            {
                for ( int i = 0; i < rank; i++ )
                    for ( int j = 0; j < HalfNN; j++ )
                        matT[j][i] = mat[i][j]; 

                diff1 = genGuesseDiffs ( matT, HalfNN, rank ); 

                if ( verbose )
                {
                    cout << "Find the first one! The practical query complexity is  " << query << " = 2^" << log2( query ) << endl;
                }

                break;
            }
        }

        //getchar();
        S1 += query;


        int S2Q = 0;
    reloop:
        for ( query = 0; query < ( 1LL << 32 ); query++ )
        {
            // the input we are waiting for
            for ( int i = 0; i < NN; i++ )
                Input2[i] = dis( gen );

            for ( int i = 0; i < HalfNN; i++ )
                V[0][i] = ( Input2[i] + Input2[ i + HalfNN ] ) % 2;

            // the difference rotations
            for ( int row = 1; row < HalfNN; row++ )
                for ( int j = row; j < HalfNN + row; j++ )
                    V[row][j % HalfNN] = V[0][ j - row ];

            for ( int i = 0; i < HalfNN; i++ )
                for ( int j = 0; j < HalfNN; j++ )
                    mat[j][i] = V[i][j];

            int rank = Gauss_F2( mat, HalfNN, HalfNN );   

            if ( rank == RR )
            {
                for ( int i = 0; i < rank; i++ )
                    for ( int j = 0; j < HalfNN; j++ )
                        matT[j][i] = mat[i][j]; 

                diff2 = genGuesseDiffs ( matT, HalfNN, rank ); 

                F2 ** MatT = new F2 * [NN];
                for ( int i = 0; i < NN; i++ )
                    MatT[i] = new F2[ NN ];

                for ( int row = 0; row < NN; row++ )
                    for ( int col = row; col < NN + row; col++ )
                        MatT[col%NN][row] = Input2[col - row];

                // [1,0,0,0] 
                // [0, 1,0,0]

                // generate the linear mapping U2
                //int RESC = 0;
                for ( int i = 0; i < HalfNN + RR; i++ )
                {
                    F2 * tempb = new F2 [ NN ];
                    for ( int j = 0; j < NN; j++ )
                        tempb[(j + i)%NN] = Input1[j];  

                    auto res = solveF2( MatT, tempb, NN, NN );

                    if ( res.size() == 0 )
                    {
                        //cout << "No solution " << endl;
                        S2Q += query;
                        //cout << "S2 " << S2 << endl;
                        //getchar();
                        //continue;
                        //return -2;
                        goto reloop;
                    }
                    else
                    {
                        ;
                        //cout << i << " " << NN << " " << res.size() << endl;
                    }

                    for ( int j = 0; j < NN; j++ )
                        U2[i][j] = res[0][j];

                    delete [] tempb;
                }

                for ( int i = 0; i < NN; i++ )
                    delete [] MatT[i];
                delete [] MatT;

                //cout << NN << " " << RESC << endl;


                //int sum = 0;
                //for ( int i = 0; i < NN; i++ )
                //    sum += Input2[i];

                if ( verbose )
                {
                    cout << "Find the second one! The practical query complexity is  " << query << " = 2^" << log2( query ) << endl;
                }
                S2Q += query;

                break;
            }
        }

        S2 += S2Q;


        //printVec( Input1, NN );
        //printVec( Input2, NN );

        /*
        int sum = 0;

        for ( int i = 0; i < NN; i++ )
        {
            cout << ( Input1[i] ^ Input2[i] );
            sum += ( Input1[i] ^ Input2[i] );
        }
        cout << endl << sum <<  endl;
        */


        //getchar();

        if ( verbose )
        {
            cout << "U2 is " << endl;
            printMatrix( U2, HalfNN + RR, NN );
        }


        //getchar();

        if ( verbose )
        {
            cout << "The first input is " << endl;
            cout << "I1 = \t "; 
            for ( int i = 0; i < NN; i++ )
                cout << static_cast<int> ( Input1[i] );
            cout << endl;

            cout << "The first output is " << endl;
            cout << "O1 = \t ";
        }
        wPRF( KEY, B, Input1, Output1 );

        if ( verbose )
        {
            cout << "The first output is " << endl;
            cout << "O1 = \t ";

            for ( int i = 0; i < TT; i++ )
                cout << static_cast<int> ( Output1[i] );
            cout << endl;
            cout << "Transform the output to the F3*-homomorphism output " << endl;
        }

        for ( int i = 0; i < TT; i++ )
            for ( int j = 0; j < NN; j++ )
                Output1[i] = ( Output1[i] + B[i][j] ) % 3;

        if ( verbose )
        {
            cout << "The F3*-homomorphism of the first output is " << endl;
            cout << "O1' = \t";
            for ( int i = 0; i < TT; i++ )
                cout << static_cast<int> ( Output1[i] );
            cout << endl;

            cout << "The second input is " << endl;
            cout << "I2 = \t "; 
            for ( int i = 0; i < NN; i++ )
                cout << static_cast<int> ( Input2[i] );
            cout << endl;

            cout << "The second output is " << endl;
            cout << "O2 = \t ";
        }

        wPRF( KEY, B, Input2, Output2 );

        if ( verbose )
        {
            for ( int i = 0; i < TT; i++ )
                cout << static_cast<int> ( Output2[i] );
            cout << endl;

            cout << "Transform the output to the F3*-homomorphism output " << endl;
        }
        for ( int i = 0; i < TT; i++ )
            for ( int j = 0; j < NN; j++ )
                Output2[i] = ( Output2[i] + B[i][j] ) % 3;

        if ( verbose )
        {
            cout << "The F3*-homomorphism of the second output is " << endl;
            cout << "O2' = \t";
            for ( int i = 0; i < TT; i++ )
                cout << static_cast<int> ( Output2[i] );
            cout << endl;

            cout << "first mid " << endl;
        }

        F2 * mid1 = new F2 [ NN ];
        MatrixMulVector2( KEY, Input1, mid1, NN );
        //printVec( mid1, NN );

        if ( verbose ) 
        {
            cout << "second mid " << endl;
        }

        F2 * mid2 = new F2 [ NN ];
        MatrixMulVector2( KEY, Input2, mid2, NN );
        //printVec( mid2, NN );


        // check the projection
        F2 * pro1 = new F2[ HalfNN + RR ];
        for ( int i = 0; i < HalfNN + RR; i++ )
            pro1[i] = 0;

        for ( int i = 0; i < HalfNN+RR; i++ )
            for ( int j = 0; j < NN; j++ )
                pro1[i] ^= U1[i][j] * mid1[j]; 

        F2 * pro2 = new F2[ HalfNN + RR ];

        for ( int i = 0; i < HalfNN + RR; i++ )
            pro2[i] = 0;

        for ( int i = 0; i < HalfNN+RR; i++ )
            for ( int j = 0; j < NN; j++ )
                pro2[i] ^= U2[i][j] * mid2[j]; 

        if ( verbose )
        {
            printVec( pro1, HalfNN + RR );
            printVec( pro2, HalfNN + RR );
        }

        delete [] pro1;
        delete [] pro2;

        //getchar();

        if ( verbose )
        {
            cout << "-----------------------------------------------------------" << endl;
            cout << "Online phase " << endl;
        }

        // guess the values in the first half, and generate the second half
        // will guess the values of first ( HalfNN - TT ) variables
        //
        F3 * tempO = new F3[TT];
        F3 ** tempB =  new F3 * [TT];
        for ( int i = 0; i < TT; i++ )
            tempB[i] = new F3[TT];

        vector< vector<F2> > solutions;

        if ( verbose )
        {
            cout << "Online phase complexity: 2 x 2^" << ( RR + HalfNN - TT ) << endl;   
        }

        //cout << diff1.size() << " " << diff2.size() << endl;
        
        vector< vector<F2>> solution1;
        vector< vector<F2>> solution2;

        for ( auto it : diff2 ) // for each difference
        //for ( int index = 0; index < ( 1 << RR ); index ++ )
        {
            // **********|vvvvvvvvvv|**********|vvvvvvvvvv
            for ( int x = 0; x < ( 1 << (HalfNN - TT) ); x++ ) // guess the values for the first HalfNN - TT 
            {
                // prepare the coefficients
                for ( int t = 0; t < TT; t++ )
                    for ( int j = 0; j < TT; j++ )
                        tempB[t][j] = ( B[t][ HalfNN - TT + j ] + ( B[t][HalfNN + HalfNN - TT + j ] * ( it[ HalfNN - TT + j] + 1 ) ) ) % 3;

                // prepare the output1
                for ( int t = 0; t < TT; t++ )
                {
                    F3 reminder = 0;

                    for ( int i = 0; i < ( HalfNN - TT ); i++ )
                    {
                        reminder = ( reminder + B[t][i] * ( ( x >> ( HalfNN - TT - 1 - i ) & 0x1 ) + 1 ) ) % 3; // + 1 to make it into F3*
                        reminder = ( reminder + B[t][i + HalfNN] * ( ( x >> ( HalfNN - TT - 1 - i ) & 0x1 ) + 1 ) * ( it[i] + 1 ) ) % 3; // F3* is multiplicative group
                    }

                    tempO[t] = mod3 ( Output1[t] - reminder ); 
                }

                auto res1 = solveF3( tempB, tempO, TT, TT ); 

                W1_before += res1.size();

                if ( res1.size() == 0 )
                    continue;

                for ( auto res : res1 )
                {
                    bool flag = true;
                    for ( int i = 0; i < TT; i++ )
                    {
                        if ( res[i] == 0 )
                        {
                            flag = false; 
                            break;
                        }
                    }
                    if ( flag == false )
                        continue;

                    vector<F2> sol ( NN );
                    for ( int i = 0; i < ( HalfNN - TT ); i++ )
                        sol[i] = x >> ( HalfNN - TT - 1 - i ) & 0x1;

                    for ( int i = 0; i < TT; i++ )
                        sol[HalfNN - TT + i] = ( res[i] - 1 );

                    for ( int i = 0; i < ( HalfNN - TT ); i++ )
                        sol[HalfNN + i] = ( ( x >> ( HalfNN - TT - 1 - i ) & 0x1 ) ^ it[i] );

                    for ( int i = 0; i < TT; i++ )
                        sol[HalfNN + HalfNN - TT + i] =  static_cast<int> ( ( res[i] - 1 ) ^ ( it[HalfNN - TT + i]) );

                    solution1.push_back( sol );
                }
            }

            // **********|vvvvvvvvvv|**********|vvvvvvvvvv
            for ( int x = 0; x < ( 1 << (HalfNN - TT) ); x++ ) // guess the values for the first HalfNN - TT 
            {
                // prepare the coefficients
                for ( int t = 0; t < TT; t++ )
                    for ( int j = 0; j < TT; j++ )
                        tempB[t][j] = ( B[t][ HalfNN - TT + j ] + ( B[t][HalfNN + HalfNN - TT + j ] * ( it[ HalfNN - TT + j] + 1 ) ) ) % 3;

                // prepare the output1
                for ( int t = 0; t < TT; t++ )
                {
                    F3 reminder = 0;

                    for ( int i = 0; i < ( HalfNN - TT ); i++ )
                    {
                        reminder = ( reminder + B[t][i] * ( ( x >> ( HalfNN - TT - 1 - i ) & 0x1 ) + 1 ) ) % 3; // + 1 to make it into F3*
                        reminder = ( reminder + B[t][i + HalfNN] * ( ( x >> ( HalfNN - TT - 1 - i ) & 0x1 ) + 1 ) * ( it[i] + 1 ) ) % 3; // F3* is multiplicative group
                    }

                    tempO[t] = mod3 ( Output2[t] - reminder ); 
                }

                auto res2 = solveF3( tempB, tempO, TT, TT ); 

                if ( res2.size() == 0 )
                    continue;

                for ( auto res : res2 )
                {
                    bool flag = true;
                    for ( int i = 0; i < TT; i++ )
                    {
                        if ( res[i] == 0 )
                        {
                            flag = false; 
                            break;
                        }
                    }

                    if ( flag == false )
                        continue;

                    vector<F2> sol ( NN );
                    for ( int i = 0; i < ( HalfNN - TT ); i++ )
                        sol[i] = x >> ( HalfNN - TT - 1 - i ) & 0x1;

                    for ( int i = 0; i < TT; i++ )
                        sol[HalfNN - TT + i] = ( res[i] - 1 );

                    for ( int i = 0; i < ( HalfNN - TT ); i++ )
                        sol[HalfNN + i] = ( ( x >> ( HalfNN - TT - 1 - i ) & 0x1 ) ^ it[i] );

                    for ( int i = 0; i < TT; i++ )
                        sol[HalfNN + HalfNN - TT + i] =  static_cast<int> ( ( res[i] - 1 ) ^ ( it[HalfNN - TT + i]) );

                    solution2.push_back( sol );
                }
            }
            // generate the linear equations
        }

        //cout << solution1.size() << " " << solution2.size() << endl;

        W1 +=  solution1.size();
        W2 +=  solution2.size();

        vector<F2> realmid1 ( mid1, mid1 + NN );
        if ( verbose )
        {
        if ( find( solution1.begin(), solution1.end(), realmid1 ) != solution1.end() )
            cout << "Find it in solution1 " << endl;
        }

        delete [] mid1;

        vector<F2> realmid2 ( mid2, mid2 + NN );

        if ( verbose )
        {
        if ( find( solution2.begin(), solution2.end(), realmid2 ) != solution2.end() )
            cout << "Find it in solution2 " << endl;
        }

        delete [] mid2;

        set< vector<F2> > Hash_table_for_w1;

        for ( auto it : solution1 )
        {
            vector<F2> s( NN, 0 );

            for ( int i = 0; i < HalfNN + RR; i++ )
                for ( int j = 0; j < NN; j++ )
                    s[i] ^= U1[i][j] * it[j];

            Hash_table_for_w1.insert( s );
        }


        if ( verbose )
        {
            cout << "After insection, the correct w' is " << endl;
        }

        F2 * corw = new F2[NN];

        int countt = 0;
        for ( auto it : solution2 )
        {
            vector<F2> s( NN, 0 );
            for ( int i = 0; i < HalfNN + RR; i++ )
                for ( int j = 0; j < NN; j++ )
                    s[i] ^= U2[i][j] * it[j];
            
            if ( Hash_table_for_w1.count( s ) > 0 )
            {
                countt += 1;

                
                //printVec( it, NN );

                for (int i = 0; i < NN; i++ )
                    corw[i] = it[i];
            }
        }

        F2 ** M2 = new F2*[NN];
        for (int i = 0; i < NN; i++ )
           M2[i] = new F2[NN]; 

        for ( int row = 0; row < NN; row++ )
            for (int col = row; col < NN + row; col++ )
                M2[row][ col % NN ] = Input2[ col - row ];

        auto res = solveF2( M2, corw, NN, NN );


        if ( verbose )
        {
            cout << "Key candidate Size " <<  res.size() << endl;
        }

        F2 * Input3 = new F2[NN];

        F3 * Output3 = new F3[TT];
        F3 * Output3Candi = new F3[TT];

        F2 * Input4 = new F2[NN];

        F3 * Output4 = new F3[TT];
        F3 * Output4Candi = new F3[TT];

        F2 * Input5 = new F2[NN];

        F3 * Output5 = new F3[TT];
        F3 * Output5Candi = new F3[TT];


        for ( int i = 0; i < NN; i++ )
            Input3[i] = dis( gen );

        for ( int i = 0; i < NN; i++ )
            Input4[i] = dis( gen );

        for ( int i = 0; i < NN; i++ )
            Input5[i] = dis( gen );

        //cout << "Input 3 " << endl;
        //printVec( Input3, NN );

        wPRF( KEY, B, Input3, Output3 );
        wPRF( KEY, B, Input4, Output4 );
        wPRF( KEY, B, Input5, Output5 );

        //printVec( Input3, NN );

        //getchar();

        F2 ** KEY1 = new F2 * [NN];
        for ( int i = 0; i < NN; i++ ) 
            KEY1[i] = new F2[NN];

        E += res.size();

        for ( auto & it : res )
        {
            for ( int i = 0; i < NN; i++ )
                for ( int j = 0; j < NN; j++ ) 
                    KEY1[i][j] = it[ ( i + j ) % NN ];

            wPRF( KEY1, B, Input3, Output3Candi );
            wPRF( KEY1, B, Input4, Output4Candi );
            wPRF( KEY1, B, Input5, Output5Candi );

            if ( EQ( Output3, Output3Candi, TT ) && 
                 EQ( Output4, Output4Candi, TT ) && 
                 EQ( Output5, Output5Candi, TT ) 
                )
            {
                //cout << "Successful " << endl;
                SUC += 1;

                if ( verbose )
                {
                    cout << "The recovered key is " << endl;
                    printVec( it, NN );
                    printVec( KEY[0], NN );
                }

                int s3 = 0, s4 = 0, s5 = 0;
                for ( int i = 0; i < NN; i++ )
                {
                    s3 += Input3[i];
                    s4 += Input4[i];
                    s5 += Input5[i];
                }

                //cout << s3 << " " << s4 << " " << s5 << endl;

                break;
            }
            //cout << count ++  << endl;
        }

        for ( int i = 0; i < HalfNN + RR; i++ )
            delete [] U1[i];
        delete [] U1;

        for ( int i = 0; i < HalfNN + RR; i++ )
            delete [] U2[i];
        delete [] U2;

        for ( int i = 0; i < NN; i++ )
            delete [] M2[i];
        delete [] M2;

        delete [] corw;

        for ( int i = 0; i < NN; i++ )
            delete [] KEY[i];
        delete [] KEY;

        for ( int i = 0; i < TT; i++ )
            delete [] B[i];
        delete [] B;

        for ( int i = 0; i < NN; i++ )
            delete [] KEY1[i];
        delete [] KEY1;

        delete [] Input1;
        delete [] Input2;
        delete [] Input3;
        delete [] Input4;
        delete [] Input5;
        delete [] Output1;
        delete [] Output2;
        delete [] Output3;
        delete [] Output4;
        delete [] Output5;
        delete [] Output3Candi;
        delete [] Output4Candi;
        delete [] Output5Candi;

        delete [] tempO;
        for ( int i = 0; i < TT; i++ )
            delete [] tempB[i];
        delete [] tempB;

        for ( int i = 0; i < HalfNN; i++ )
            delete [] mat[i];
        delete [] mat;

        for ( int i = 0; i < RR; i++ )
            delete [] matT[i];
        delete [] matT;

        for ( int i = 0; i < HalfNN; i++ )
            delete [] V[i];
        delete [] V;

        if ( COUNT == TESTNUM )
            break;

        const int barWidth = 50; 

         float progress = static_cast<float>(COUNT) / TESTNUM;
            int pos = barWidth * progress;

        std::string bar;
        bar.reserve(barWidth + 10); 
        bar = "[";
        for (int p = 0; p < barWidth; ++p) {
            bar += (p < pos) ? '=' : ' ';
        }
        bar += "] " + std::to_string(static_cast<int>(progress * 100)) + "%";

        std::cout << "\r" << bar << std::flush;
    }

    cout << "Done!" << endl;

    cout << endl;

    cout << "The average complexity over " << TESTNUM << " tests: " << endl;

    cout << "Sampling1  " << log2( S1 / TESTNUM ) << endl;
    cout << "Sampling2  " << log2( S2 / TESTNUM ) << endl;
    cout << "W1_before  " << log2( W1_before / TESTNUM ) << endl;
    cout << "W1         " << log2( W1 / TESTNUM ) << endl;
    cout << "W2         " << log2( W2 / TESTNUM ) << endl;
    cout << "Exhaustive " << log2( E / TESTNUM ) << endl;
    cout << "Success    " << SUC << endl;

    //for ( int i = 0; i < 100; i++ )
    //    cout << i << " " << Value[i] << endl;

}

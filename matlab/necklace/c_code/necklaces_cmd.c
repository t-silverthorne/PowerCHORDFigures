//------------------------------------------------------
// NECKLACES, LYNDON WORDS, and RELATIVES
// Programmed by: Joe Sawada 2010-2014
//------------------------------------------------------
#include <stdio.h>
#define MAX 99
#define TRUE 1
#define FALSE 0

typedef struct cell {
    int next,prev;
} cell;

typedef struct element {
    int s, v;
} element;

//-------------------------------------------------------------
// GLOBAL VARIABLES
//-------------------------------------------------------------
int N,K,D,M,type,head,total=0,NECK=0, LYN=0;
int UNRESTRICTED=0,DENSITY=0,CONTENT=0,FORBIDDEN=0,BRACELET=0,UNLABELED=0,CHORD=0,LIE=0,CHARM=0,DB=0;
int a[MAX], p[MAX], b[MAX], f[MAX], fail[MAX], num[MAX], run[MAX], num_map[MAX], charm[2*MAX+2];
int pos[MAX][MAX], split[MAX][MAX], d[MAX][MAX], match[MAX][MAX];
cell avail[MAX];
int nb = 0; // number of blocks
element B[MAX]; // run length encoding data structure
char PRIME[MAX]; // relatively prime array for charm bracelets

//-------------------------------------------------------------
void Input(int argc, char *argv[]) {
    int i, j, n_digit;
    
    sscanf(argv[1], "%d", &type);
    
    if (type == 1 || type == 21) UNRESTRICTED = 1;		// SLOANE A000031 and A001037
    if (type == 2 || type == 22) DENSITY = 1;
    if (type == 3 || type == 6 || type == 9 || type == 23 || type == 26) CONTENT = 1;
    if (type == 4 || type == 24) FORBIDDEN = 1;
    if (type == 5 || type == 6 || type == 8 || type == 9 ||
        type == 25 || type == 26) BRACELET = 1;			// SLOANE A000029 and A001371
    if (type == 7 || type == 27) UNLABELED = 1;			// SLOANE A000013 and A000048
    if (type == 8 || type == 9) CHARM = 1;				// SLOANE A002729
    if (type == 10) CHORD = 1;							// SLOANE A007769
    if (type == 11) LIE = 1;
    if (type == 12) UNRESTRICTED = DB = 1;
    
    if (type <= 20) NECK=1;
    else LYN=1;
    
    //---------------------------------------
    if (CHORD) printf("\n  ENTER number of chords n: ");
    sscanf(argv[2], "%d", &N);
    //---------------------------------------
    K=2;
    if (!UNLABELED && !CHORD) {
        sscanf(argv[3], "%d", &K);
    }
    //---------------------------------------
    if (DENSITY) {
        sscanf(argv[4], "%d", &D);
    }
    //---------------------------------------
    if (CONTENT) {
        i = 1;
        for (j=1; j<=K; j++) {
            scanf("%d", &n_digit);
            
            // The fixed content algo expects each num[i] > 0, but to
            // allow 0 entries, we create a number map for printing
            if (n_digit != 0) {
                num[i] = n_digit;
                num_map[i] = j;
                i++;
            }
        }
        K = i-1;
    }
    //---------------------------------------
    if (FORBIDDEN) {
        printf("  ENTER length of forbidden sequence: "); scanf("%d", &M);
        printf("\n");
        for (i=1; i<=M; i++) {
            printf("  ENTER (character < k) forb[%d]: ",i);
            scanf("%d", &f[i]);
        }
    }
    printf("\n");
}
//-------------------------------------------------------------
//  SUBROUTINES FOR CHARM BRACELETS
//-------------------------------------------------------------
// finds the index of the minimal necklace rotation
int MinNecklace(){
    int i,k,p;
    
    // create necklace of string b[n]
    k = 0; charm[2*N+1] = -1;
    while (k < 2*N) {
        i = k +2;
        p = 1;
        while (charm[i-p] <= charm[i]) {
            if (charm[i-p] < charm[i]) p = i-k;
            i++;
        }
        if (p * ((i-k-1)/p) >=N) break;
        do {
            k = k+p;
        } while ( k < i-p);
    }
    
    return k+1;
}
//-------------------------------------------------------------
int Gcd(int x, int y){
    int t;
    
    while( y != 0 ) {
        t = y;
        y = x % y;
        x = t;
    }
    return x;
}
//-------------------------------------------------------------
int IsCharm(){
    int i,j,offset;
    
    // skip N-1 since it's taken care in the bracelet generation algorithm
    for(i=2; i<N-1; i++){
        
        // only consider numbers relatively prime to N
        if ( PRIME[i] ) {
            for(j=0; j<N; j++) charm[j*i % N + 1] = charm[j*i % N + N+1] = a[j+1];
            offset = MinNecklace();
            
            for (j=1; j<=N; j++){
                if (a[j] < charm[offset + j-1]) break;
                else if (a[j] > charm[offset + j-1]) return FALSE;
            }
        }
    }
    return TRUE;
}
//-------------------------------------------------------------
//-------------------------------------------------------------
void Print() {
    int j;
    
    if (CHORD) {
        for (j=0; j<2*N; j++) printf("%d ", a[j]);
        total++;
        printf("\n");
    }
    else {
        if (!CHARM || (CHARM && IsCharm())) {
            for (j=1; j<=N; j++) {
                if (CONTENT) printf("%d", num_map[a[j]]-1);
                else printf("%d", a[j]);
            }
            total++;
            printf("\n");
        }
    }
}
//-------------------------------------------------------------
void PrintDB(int p) {
    int j;
    
    for (j=1; j<=p; j++) printf("%d", a[j]);
    total+=p;
}
//-------------------------------------------------------------
void PrintD(int p) {
    int i,j,next,end,min;
    
    /* Determine minimum position for next bit */
    next =  (D/p)*a[p] + a[D%p];
    if (next < N) return;
    
    /* Determine last bit */
    min = 1;
    if ((next == N) && (D%p != 0)) {
        min =  b[D%p]+1;
        p = D;
    }
    else if ((next == N) && (D%p == 0)) min = b[p];
    end = N;
    for( b[D]=min; b[D]<K; b[D]++ ) {
        i = 1;
        if ( LYN && (N%a[p] == 0) && (a[p] != N)) {}
        else {
            for(j=1; j<=end; j++) {
                if (a[i] == j) {
                    printf("%d",b[i]);
                    i++;
                }
                else printf("0");
            }
            printf("\n");
            total++;
        }
        p = D;
    }
}
//-------------------------------------------------------------
void PrintBracket(int start, int end) {
    
    if (start == end)  printf("%d", a[start]);
    else {
        printf("[ ");
        PrintBracket(start,split[start][end]-1);
        printf(" , ");
        PrintBracket(split[start][end],end);
        printf(" ]");
    }
}
/*------------------------------------------------------------*/
// SUBROUTINES
/*------------------------------------------------------------*/
void SetMatch(int t) {
    int j;
    
    for (j=t; j<M; j++) {
        if (match[j][t-1] && f[M-j+t] == a[t]) match[j][t] = TRUE;
        else match[j][t] = FALSE;
    }
}
/*------------------------------------------------------------*/
int CheckSuffix(int s) {
    
    while(s > 0) {
        if ( match[M-s][M-s] == TRUE ) return FALSE;
        else s = fail[s];
    }
    return TRUE;
}
/*------------------------------------------------------------*/
void Remove(int i) {
    int p,n;
    
    if (i == head) head = avail[i].next;
    p = avail[i].prev;
    n = avail[i].next;
    avail[p].next = n;
    avail[n].prev = p;
}
/*------------------------------------------------------------*/
void Add(int i) {
    int p,n;
    
    p = avail[i].prev;
    n = avail[i].next;
    avail[n].prev = i;
    avail[p].next = i;
    if (avail[i].prev == K+1) head = i;
}

/*------------------------------------------------------------*/
void InitListC() {
    int i;
    
    for (i=1; i<N*2; i++) {
        avail[i].next = i+1;
        avail[i].prev = i-1;
    }
    head = 1;
}
/*------------------------------------------------------------*/
void RemoveC( int i) {
    int p,n;
    
    if (avail[i].prev == 0) {
        head = avail[i].next;
        avail[head].prev = 0;
    }
    else  {
        p = avail[i].prev;
        n = avail[i].next;
        avail[p].next = n;
        avail[n].prev = p;
    }
    avail[i].next = i;
}
/*------------------------------------------------------------*/
void AddC( int i) {
    int p,n;
    
    if (avail[i].prev == 0) {
        avail[head].prev = i;
        avail[i].next = head;
        head = i;
    }
    else {
        p = avail[i].prev;
        n = avail[p].next;
        avail[i].next = n;
        avail[n].prev = i;
        avail[p].next = i;
    }
}
/*------------------------------------------------------------*/
// return -1 if reverse smaller
// return 0 if reverse larger
// return 1 if equal
/*------------------------------------------------------------*/
int CheckRev( int t, int i ) {
    int j;
    
    for (j=i+1; j<=(t+1)/2; j++) {
        if      (a[j] < a[t-j+1]) return(0);
        else if (a[j] > a[t-j+1]) return(-1);
    }
    return(1);
}
/*------------------------------------------------------------*/
// SUBROUTINES FOR FIXED CONTENT (CHARM) BRACELTS
/*------------------------------------------------------------*/
// return -1 if reverse smaller
// return 0  if equal
// return 1  if reverse is larger
/*------------------------------------------------------------*/
int CheckRevF() {
    int j;
    
    j = 1;
    while (B[j].v == B[nb-j+1].v && B[j].s == B[nb-j+1].s && j<= nb/2) j++;
    if (j > nb/2) return 0;
    if (B[j].s < B[nb-j+1].s) return 1;
    if (B[j].s > B[nb-j+1].s) return -1;
    if (B[j].v < B[nb-j+1].v && B[j+1].s < B[nb-j+1].s) return 1;
    if (B[j].v > B[nb-j+1].v && B[j].s < B[nb-j].s) return 1;
    return -1;
}
/*------------------------------------------------------------*/
void UpdateRunLength(int v) {
    
    if (B[nb].s == v) B[nb].v = B[nb].v + 1;
    else {
        nb++;
        B[nb].v = 1;
        B[nb].s = v;
    }
}
/*------------------------------------------------------------*/
void RestoreRunLength() {
    
    if (B[nb].v == 1) nb--;
    else B[nb].v = B[nb].v - 1;
}
/*-----------------------------------------------------------*/
// UNRESTRICTED
/*-----------------------------------------------------------*/
void Gen(int t, int p) {
    int j;
    
    if(t>N) {
        if (DB && N%p==0) PrintDB(p);
        else if ((NECK && N%p==0) || (LYN  && N==p))  Print();
    }
    else {
        a[t] = a[t-p];
        Gen(t+1,p);
        for(j=a[t-p]+1; j<=K-1; j++) {
            a[t] = j;
            Gen(t+1,t);
        }
    }
}
/*-----------------------------------------------------------*/
// FIXED DENSITY
/*-----------------------------------------------------------*/
void GenD(int t,int p) {
    int i,j,max,tail;
    
    if (t >= D-1) PrintD(p);
    else {
        tail = N - (D - t) + 1;
        max = a[t-p+1] + a[p];
        if (max <=tail) {
            a[t+1] = max;
            b[t+1] = b[t-p+1];
            
            GenD(t+1,p);
            for (i=b[t+1] +1; i<K; i++) {
                b[t+1] = i;
                GenD(t+1,t+1);
            }
            tail = max-1;
        }
        for(j=tail; j>=a[t]+1; j--) {
            a[t+1] =  j;
            for (i=1; i<K; i++) {
                b[t+1] =  i;
                GenD(t+1,t+1);
            }
        }
    }
}
/*-----------------------------------------------------------*/
// FIXED CONTENT
/*-----------------------------------------------------------*/
void GenA(int t, int p, int s) {
    int j, s2;
    
    if (num[K] == N-t+1) {
        if (NECK && (num[K] == run[t-p]) && (N%p == 0)) Print();
        else if (LYN && (num[K] == run[t-p]) && (N==p)) Print();
        else if  (num[K] > run[t-p]) Print();
    }
    else if (num[1] != N-t+1) {
        j = head;
        s2 = s;
        while( j >= a[t-p]) {
            
            run[s] = t-s;
            a[t] = j;
            
            num[j]--;
            if (num[j] == 0) Remove(j);
            
            if (j != K)  s2 = t+1;
            if (j == a[t-p]) GenA(t+1,p,s2);
            else GenA(t+1,t,s2);
            
            if (num[j] == 0) Add(j);
            num[j]++;
            
            j = avail[j].next;
        }
        a[t] = K;
    }
}
/*------------------------------------------------------------*/
// FORBIDDEN SUBSTRINGS
/*------------------------------------------------------------*/
void GenF(int t, int p, int s) {
    int j,q;
    
    if (t > N) {
        if (NECK && N%p == 0 && CheckSuffix(s))  Print();
        if (LYN  && N==p     && CheckSuffix(s))  Print();
    }
    else {
        a[t] = a[t-p];
        q = d[s][a[t]];
        if (t<M) SetMatch(t);
        if (q!=M) GenF(t+1,p,q);
        
        for(j=a[t-p]+1; j<=K-1; j++) {
            a[t] = j;
            q = d[s][j];
            if (t < M) SetMatch(t);
            if (q != M) GenF(t+1,t,q);
        }
    }
}
/*-----------------------------------------------------------*/
// BRACELETS
/*-----------------------------------------------------------*/
void GenB( int t, int p, int r, int u, int v, int RS ) {
    int j, rev;
    
    if (t-1 > (N-r)/2 + r) {
        if (a[t-1] > a[N-t+2+r]) RS = FALSE;
        else if (a[t-1] < a[N-t+2+r]) RS = TRUE;
    }
    if (t > N)  {
        if ((RS == FALSE) && ((NECK && N%p == 0) || (LYN && N==p))) Print();
    }
   	else {
        a[t] = a[t-p];
        if (a[t] == a[1]) v++;
        else v = 0;
        if ((u == -1) && (a[t-1] != a[1])) u = r =  t-2;
        
        if ((u != -1) && (t == N) && (a[N] == a[1])) {}
        else if (u == v) {
            rev = CheckRev(t, u);
            if (rev == 0) GenB(t+1,p,r,u,v,RS);
            if (rev == 1) GenB(t+1,p,t,u,v,FALSE);
        }
        else GenB(t+1,p,r,u,v,RS);
        for (j = a[t-p]+1; j <= K-1; ++j) {
            a[t] = j;
            GenB(t+1,t,r,u,0,RS);
        }
   	}
}
/*-----------------------------------------------------------*/
// FIXED CONTENT BRACELETS
/*-----------------------------------------------------------*/
void GenBF(int t, int p, int r, int z, int b, int RS) {
    int j,z2,p2,c;
    int i,count,order[20];  // used for lex generation
    
    // Incremental comparison of a[r+1...n] with its reversal
    if (t-1 > (N-r)/2 + r) {
        if (a[t-1] > a[N-t+2+r]) RS = FALSE;
        else if (a[t-1] < a[N-t+2+r]) RS = TRUE;
    }
    
    // Termination condition - only characters k remain to be appended
    if (num[K] == N-t+1) {
        if (num[K] > run[t-p]) p = N;
        if (num[K] > 0 && t != r+1 && B[b+1].s == K && B[b+1].v > num[K]) RS = TRUE;
        if (num[K] > 0 && t != r+1 && (B[b+1].s != K || B[b+1].v < num[K])) RS = FALSE;
        if ((RS == FALSE) && ((NECK && N%p == 0) || (LYN && N==p))) Print();
    }
    // Recursively extend the prenecklace - unless only 0s remain to be appended
    else if (num[1] != N-t+1) {
        j = head;
        
        // generate in lex order
        i = 0;
       	while (j >= a[t-p]) {
            order[++i] = j;
            j = avail[j].next;
       	}
        
        //while( j >= a[t-p]) {
        
        for (count = i; count>=1; count--) {
            j = order[count];
            
            run[z] = t-z;
            UpdateRunLength(j);
            num[j]--;
            if (num[j] == 0) Remove(j);
            a[t] = j;
            z2 = z;
            if (j != K) z2 = t+1;
            p2 = p;
            if (j != a[t-p]) p2 = t;
            c = CheckRevF();
            if (c == 0) GenBF(t+1,p2,t,z2,nb,FALSE);
            if (c == 1) GenBF(t+1,p2,r,z2,b,RS);
            if (num[j] == 0) Add(j);
            num[j]++;
            RestoreRunLength();
            //j = avail[j].next;
        }
        a[t] = K;
    }
}
/*-----------------------------------------------------------*/
// UNLABELED
/*-----------------------------------------------------------*/
void GenU(int t, int p, int c) {
    
    if (t>N) {
        if ((NECK && N%p==0) || (LYN  && N==p))  Print();
    }
    else {
        if (a[t-c] == 0) {
            if (a[t-p] == 0) {
                a[t] =0; GenU(t+1,p,t);
            }
            a[t] = 1;
            if (a[t-p] == 1) GenU(t+1,p,c);
            else GenU(t+1,t,c);
        }
        else {
            a[t] = 0; GenU(t+1,p,c);
        }
    }
}
/*-----------------------------------------------------------*/
// CHORD DIAGRAMS
/*-----------------------------------------------------------*/
void GenRestC(int s, int e, int v) {
    
    if (s == N*2) Print();
    else if (e < N*2) {
        if (e-s >= v && N*2-e+s >= v) {
            a[s] = e - s; a[e] = N*2 - a[s];
            RemoveC(s); RemoveC(e);
            GenRestC(head,avail[head].next,v);
            AddC(e); AddC(s);
        }
        if (N*2-avail[e].next+s >= v) GenRestC(s,avail[e].next,v);
    }
}
/*-----------------------------------------------------------*/
void GenC2(int s, int t, int p, int T, int P, int len, int num_sec) {
    int e, next, min;
    
    min  = pos[len][t-p];
    if (len == N && s > N) {
        if (T == N) Print();
    }
    else if (t == 1 && s > P) {
        if (len < N) GenC2(head,1,1,T,P,len+1,0);
    }
    else if (s > (num_sec+1)*P ) {
        pos[len][t] = P;
        if (min == P)  GenC2(s,t+1,p,T,P,len,num_sec+1);
        else  GenC2(s,t+1,t,T,P,len,num_sec+1);
    }
    else if (s == N*2) {
        if (T == N) Print();
        else if (min != P) GenRestC(head,avail[head].next,len+1);
        else if (t%p == 0 && len < N) GenC2(head,1,1,T,N*2*p/t,len+1,0);
    }
    else {
        e = (s + len) % (N*2);
        if (s % P >= min && avail[e].next != e) {
            next = avail[s].next;
            if (next == e) next = avail[e].next;
            a[s] = (N*2+ e -s) % (N*2);
            a[e] = N*2 - a[s];
            RemoveC(s); RemoveC(e);
            pos[len][t] = s % P;
            if (s % P == min) GenC2(next,t+1,p,T+1,P,len,num_sec);
            else GenC2(next,t+1,t,T+1,P,len,num_sec);
            AddC(e); AddC(s);
        }
        GenC2(avail[s].next,t,p,T,P,len,num_sec);
    }
}
/*-----------------------------------------------------------*/
void GenC(int s, int t, int p, int v, int last) {
    int next,e,min,p2;
    
    min = last + pos[v][t-p];
    if (min == N*2 && t%p == 0) {
        if (t == N) Print();
        else GenC2(head,1,1,t,N*2*p/t,v+1,0);
    }
    else if (min < N*2 && s < N*2) {
        e = (s + v) % (N*2);
        if (s >= min && avail[e].next != e) {
            
            next = avail[s].next;
            if (next == e) next = avail[e].next;
            a[s] = v; a[e] = N*2 - v;
            RemoveC(s); RemoveC(e);
            pos[v][t] = s-last;
            
            p2 = p;
            if (s != min) p2 = t;
            GenC(next,t+1,p2,v,s);
            if (s + pos[v][t+1-p2] < N*2)  GenRestC(head, avail[head].next,v+1);
            AddC(e); AddC(s);
        }
        GenC(avail[s].next,t,p,v,last);
    }
}
/*-----------------------------------------------------------*/
// LIE BRACKETS (BASIS)
/*-----------------------------------------------------------*/
void GenL(int t) {
    int i,j,q[50];
    
    if(t>N) {
        if (N == p[1]) { PrintBracket(1,N); printf("\n");  }
    }
    else {
        for(i=1; i<=N; i++) q[i] = p[i];
        for(j=a[t-p[1]]; j<=K-1; j++) {
            a[t] = j;
            for (i=1; i<=t-1; i++) {
                if (a[t] < a[ t-p[i] ]) p[i] = 0;
                if (a[t] > a[ t-p[i] ]) p[i] = t-i+1;
            }
            for (i=t-1; i >=1; i--) {
                if (p[i+1] == t-i) split[i][t] = i+1;
                else split[i][t] = split[i+1][t];
            }
            GenL(t+1);
            for(i=1; i<=N; i++) p[i] = q[i];
        }
    }
}
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
void Init() {
    int i,j;	
    
    if (CHARM) for (i=1; i<N; i++) PRIME[i] = (Gcd(N, i) == 1);
    
    a[0] = a[1] = 0;
    //----------------
    if (UNRESTRICTED)  Gen(1,1);
    //----------------
    if (DENSITY) {				
        for(j=0; j<=D; j++) a[j] = 0;		
        if (D == 0) {
            if (NECK) {
                for (j=1; j<=N; j++) printf("0");
                printf("\n");
                total = 1; 
            }
        } 
        else if (D == 1) {
            for (i=1; i<K; i++) {
                for (j=1; j<N; j++) printf("0");
                printf("%d \n",i);
            }
            total = K-1; 
        } 
        else {			
            a[D] = N; 	
            for(j=N-D+1; j>=(N-1)/D + 1; j--) {
                a[1] = j;
                for (i=1; i<K; i++) {
                    b[1] = i;
                    GenD(1,1);		
                }
            } 
        }
    }
    //----------------
    if (CONTENT) {
        for (j=K+1; j>=0; j--) {
            avail[j].next = j-1;
            avail[j].prev = j+1;
        }
        head = K;
        
        for (j=1; j<=N; j++) {
            a[j] = K;
            run[j] = 0;
        }		
        a[1] = 1;
        num[1]--;
        if (num[1] == 0) Remove(1);
        
        if (BRACELET) {
            B[0].s = 0;
            UpdateRunLength(1);
            GenBF(2,1,1,2,1,FALSE);
        }
        else GenA(2,1,2);
    }
    //----------------
    if (FORBIDDEN) {
        /* Pre process */
        for (j=1; j<M; j++) match[j][0] = TRUE;
        
        /* Failure Function */
        fail[1] = 0;
        for (j=2; j<=M; j++) {
            i = fail[j-1];
            while (f[j] != f[i+1] && i > 0 ) i = fail[i];
            if (f[j] != f[i+1] && i == 0) fail[j] = 0;
            else fail[j] = i+1;
        }	
        
        /* Transition Function */
        for (j=1; j<=M; j++) d[j-1][f[j]] = j;
        for (j=0; j<K; j++) if (j != f[1]) d[0][j] = 0;
        for (j=1; j<=M; j++) 
            for(i=0; i<K; i++) 
                if (i != f[j+1]) d[j][i] = d[ fail[j] ][i]; 
        
        GenF(1,1,0);
    }
    //----------------
    if (UNLABELED) GenU(2,1,1); 
    //----------------
    if (BRACELET && !CONTENT) GenB(1,1,1,-1,0,FALSE);
    //----------------
    if (CHORD) {
        InitListC();
        for (i=1; i<N; i++) pos[i][0] = 0;	
        for (i=1; i<N; i++) {
            a[0] = i;
            a[i] = N*2-i;
            RemoveC(i);
            GenC(head,1,1,i,0);   
            GenRestC(head, avail[head].next, i+1);
            AddC(i);
        }
        for(i=0; i<N*2; i++) a[i] = N; 
        Print();
    }
    //----------------
    if (LIE) {
        p[0] = 0;
        for (i=1; i<=N; i++)  p[i] = 1; 
        GenL(1);
    }
}
//--------------------------------------------------------------------------------
int main(int argc, char *argv[]) {
    Input(argc, argv);
    Init();
    if (!LIE) printf("\nTotal = %d\n\n", total);
    return 0;
}
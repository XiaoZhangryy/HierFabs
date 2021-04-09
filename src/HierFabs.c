#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>

// Loss function
double Loss(double *y, double *xb, int model, int *param, int *stautus)
{
    int i;
    int n = param[0];
    double temp, val = 0.0;

    if(model == 1) {
        for (i = 0; i < n; ++i){
            temp = y[i] - xb[i];
            val += temp*temp;
        }
        val /= 2.0*n;
    } else if (model == 2) {
        // for cox model, y is ascendant sorted
        temp = 0.0;
        for (i = n-1; i >= 0; --i){
            temp += exp(xb[i]);
            if (stautus[i]) val -= xb[i] - log(temp);
        }
    }

    return val;
}

double calculate_bic(double *score, int *param, int model, double gamma)
{
    int n = param[0];
    int q = param[2];
    int df = param[3];

    if (model == 1) {
        return 2*n*(*score) + df*log(n) + 2*gamma*lchoose(q, df);
    } else if (model == 2) {
        return 2*(*score) + df*log(n) + 2*gamma*lchoose(q, df);
    } else {
        return 2*n*(*score) + df*log(n) + 2*gamma*lchoose(q, df);
    }
}

// Heap sort for an double array, from big to small 
void swap(double *arr, int *set,int i,int j){
    double temp;
    temp = arr[i];
    arr[i] = arr[j];
    arr[j] = temp;
    int t;
    t = set[i];
    set[i] = set[j];
    set[j] = t;
}

void heapify(double *arr, int *set, int i, int n){
    int j=2*i+1;
    double t=arr[i];
    int temp=set[i];
    while(j<=n){
        if(j<n && arr[j]<arr[j+1]) j++;
        if(t<arr[j]){
            arr[i]=arr[j];
            set[i]=set[j];
            i=j;
            j=2*i+1;
        }
        else break;
    }
    arr[i]=t;
    set[i]=temp;
}

void HeapSortDouble(double *arr, int *set, int n){
    int i;
    for(i=n/2-1; i>=0; i--)// build heap
        heapify(arr, set, i, n-1);
    for(i=n-1; i>=1; i--){
        swap(arr, set, 0,i);
        heapify(arr, set, 0, i-1);
    }
}

// backward qualify for weak hierarchy, judge if a main effect k can vanish 
int b_qlf(int k, int *param, int *active, int *parent, int *child)
{
    int i, j, p1, p2;
    int ns = param[3], nm = param[4], c = child[k];

    for (i = nm; i < ns; ++i)
    {
        p1 = parent[2*i];
        p2 = parent[2*i+1];
        if ( (p1 == k) || (p2 == k) )
        {
            p1 = p1+p2-k;
            c--;
            for (j = 0; j < nm; ++j)
            {
                if (active[j] < p1) continue;
                if (active[j] > p1) return 0;
                break;
            }
        }
        if (c == 0) break;
    }
    return 1;
}

void Pre_work(double *x, double *y, double *weight, int *param, 
    double *meanx, double *sdx, int isq)
{
    // y has been centralized, while x hasn't.
    int n = param[0];
    int p = param[1];
    int i, j, k, c;
    double *x_i, *x_j, temp, var, cor;

    // mean and sd of X
    for (c = 0; c < p; ++c)
    {
        cor = 0.0;
        var = 0.0;
        x_i = x + c*n;
        for (k = 0; k < n; ++k)
        {
            cor += x_i[k];
            var += x_i[k]*x_i[k];
        }
        meanx[c] = cor/n;
        sdx[c] = sqrt( (var - n*meanx[c]*meanx[c])/(n-1) );
    }
    // isq = 1 if quadratic term is considered, and isq = 0 while not.
    for (i = 0; i < p; ++i)
    {
        x_i = x + i*n;
        for (j = i+1-isq; j < p; ++j)
        {
            x_j = x + j*n;
            cor = 0.0;
            var = 0.0;
            for (k = 0; k < n; ++k) 
            {
                temp = x_i[k]*x_j[k];
                cor += temp;
                var += temp*temp;
            }
            meanx[c] = cor/n;
            sdx[c] = sqrt( (var - n*meanx[c]*meanx[c])/(n-1) );
            c++;
        }
    }
}

// update sparse index
void usi(int *s_i, int *s_j, int *tt_b, int *tt_a, int *act, int ns, int iter)
{
    int i;
    s_i += *tt_a;
    s_j += *tt_a;
    for (i = 0; i < ns; ++i)
    {
        s_i[i] = act[i];
        s_j[i] = iter;
    }
    *tt_b = *tt_a;
    *tt_a += ns;
}

// unmormalized_beta
void unmormalized(double *beta, double *unbeta, int *parent, int *indexi, 
    int *indexj, int iter, double *intercept, double *meanx, double *sdx, 
    double meany, int total_active, int model)
{
    int i, ell;

    // Eliminate standard variance
    for (ell = 0; ell < total_active; ++ell)
        unbeta[ell] = beta[ell] / sdx[ indexi[ell] ];

    // calculate intercept.
    i = 0;
    ell = 0;
    while( i < iter )
    {
        if (model == 1) {
            intercept[i] = meany;
        } else if (model == 2) {
            intercept[i] = 0.0;
        } else {
            intercept[i] = meany;
        }
        while( *indexj == i )
        {
            intercept[i] -= unbeta[ell]*meanx[ *indexi ];
            ell++;
            indexi++;
            indexj++;
        }
        i++;
    }
}

void der(double *x, double *y, double *d, double *meanx, double *sdx, int model, 
    int *param, double *residual, double *xb, int isq, int *stautus, int *active,
    int hierarchy) 
{
    // calculate the 1st order Taylor Formula.
    int i, j, k, c, m, ell=0;
    int n = param[0], p = param[1], nm = param[4];
    double *x_i, *x_j, temp;

    if(model == 1){
        temp = 0.0;
        for (i = 0; i < n; ++i) {
            residual[i] = y[i] - xb[i];
            temp += residual[i];
        }
        for (c = 0; c < p; ++c) {
            x_i = x + c*n;
            d[c] = 0.0;
            for (j = 0; j < n; ++j)
                d[c] -= x_i[j] * residual[j];
            d[c] = (d[c] + meanx[c]*temp)/(n*sdx[c]);
        }
        // isq = 1 if quadratic terms are considered, and isq = 0 while not.
        if (hierarchy)
        {
            // Strong hierarchy
            if (!nm) return;
            for (i = 0; i <= active[nm-1]; ++i){
                x_i = x + i*n;
                if (i < active[ell])
                {
                    m = (2*(p+isq)-i)*(i+1)/2-i-1;
                    for (j = ell; j < nm; ++j) {
                        x_j = x + active[j]*n;
                        c = m + active[j];
                        d[c] = 0.0;
                        for (k = 0; k < n; ++k)
                            d[c] -= x_i[k]*x_j[k]*residual[k];
                        d[c] = (d[c] + meanx[c]*temp)/(n*sdx[c]);
                    }
                } else {
                    j = i-isq+1;
                    c = (2*(p+isq)-i)*(i+1)/2-isq;
                    for (; j < p; ++j) {
                        x_j = x + j*n;
                        d[c] = 0.0;
                        for (k = 0; k < n; ++k)
                            d[c] -= x_i[k]*x_j[k]*residual[k];
                        d[c] = (d[c] + meanx[c]*temp)/(n*sdx[c]);
                        c++;
                    }
                    ell++;
                }
            }
        } else {
            // Weak hierarchy
            if (!nm) return;
            for (i = 0; i < p; ++i){
                x_i = x + i*n;
                for (j = i-isq+1; j < p; ++j) {
                    x_j = x + j*n;
                    d[c] = 0.0;
                    for (k = 0; k < n; ++k)
                        d[c] -= x_i[k]*x_j[k]*residual[k];
                    d[c] = (d[c] + meanx[c]*temp)/(n*sdx[c]);
                    c++;
                }
            }
        }
    } else if(model == 2) {
        double theta = 0.0;
        int r;

        for (i = 0; i < param[2]; ++i) d[i] = 0.0;

        for (i = n-1; i >= 0; --i) {
            theta += exp(xb[i]);
            if (stautus[i]) {
                for (j = i; j < n; ++j) {
                    x_i = x;
                    temp = exp(xb[j])/theta;
                    for (c = 0; c < p; ++c) {
                        d[c] -= temp * (x_i[i]-x_i[j]) / sdx[c];
                        x_i += n;
                    }

                    if (hierarchy) {
                        m=0;
                        for (k = 0; k <= active[nm-1]; ++k){
                            x_i = x + k*n;
                            if (k < active[m])
                            {
                                r = (2*(p+isq-1)-k)*(k+1)/2;
                                for (ell = m; ell < nm; ++ell)
                                {
                                    c = r + active[ell];
                                    x_j = x + active[ell]*n;
                                    d[c] -= temp*(x_i[i]*x_j[i]-x_i[j]*x_j[j]) / sdx[c];
                                }
                            } else {
                                ell = k-isq+1;
                                x_j = x + ell*n;
                                c = (2*(p+isq)-k)*(k+1)/2-isq;
                                for (; ell < p; ++ell, x_j += n, ++c) {
                                    d[c] -= temp*(x_i[i]*x_j[i]-x_i[j]*x_j[j]) / sdx[c];
                                }
                                m++;
                            }
                        }
                    } else {
                        for (k = 0; k < p; ++k){
                            x_i = x + k*n;
                            ell = k-isq+1;
                            x_j = x + ell*n;
                            for (; ell < p; ++ell) {
                                d[c] -= temp*(x_i[i]*x_j[i]-x_i[j]*x_j[j]) / sdx[c];
                                c++;
                                x_j += n;
                            }
                        }
                    }
                }
            }
        }
    }
}

// take an initial step
void Initial(x, y, xb, beta, d, w, param, model, eps, lambda, losses, direction, 
    active, bic, df, meanx, sdx, residual, isq, stautus, hierarchy, gamma)
double *x, *y, *xb, *beta, *d, *w, eps, *lambda, *losses, *bic;
double *meanx, *sdx, *residual, gamma;
int *param, model, *direction, *active, *df, isq, *stautus, hierarchy;
{
    // d doesn't need initial value
    int i, k=0;
    int n = param[0], p = param[1];
    double value, temp, *xi;

    // calculate the derivative
    der(x, y, d, meanx, sdx, model, param, residual, xb, isq, stautus, active,
        hierarchy);

    // find a forward direction
    temp = 0.0;
    for (i = 0; i < p; ++i){
        value = fabs(d[i])/w[i];
        if (value > temp){
            temp = value;
            k = i;
        }
    }

    // calculate increment, update xb, lambda, beta, losses, direction, bic, df
    // param[3] (nonzero elements for step t), param[4] (nonzero main effect)
    value = eps/w[k];
    if (d[k] > 0.0) value *= -1.0;
    *beta = value;

    temp = Loss(y, xb, model, param, stautus);
    xi = x + k*n;
    value /= sdx[k];
    for (i = 0; i < n; ++i)
        xb[i] += (xi[i]-meanx[k])*value;
    *losses = Loss(y, xb, model, param, stautus);
    *lambda = (temp - *losses)/eps;

    *direction = 1;
    param[3]   = 1;
    param[4]   = 1;
    *active    = k;
    *bic       = calculate_bic(losses, param, model, gamma);
    *df        = 1;
}

// try a backward step
int Backward(x, y, xb, beta, d, betanew, weight, param, model, eps, lambda, 
    losses, parent, xi, back, hierarchy, active, direction, bic, df, 
    child, meanx, sdx, residual, isq, stautus, gamma)
double *x, *y, *xb, *beta, *d, *betanew, *weight, eps, *lambda, *losses, xi;
double *bic, *meanx, *sdx, *residual, gamma;
int *param, model, *parent, back, hierarchy, *active, *direction;
int *df, *child, isq, *stautus;
{
    if (!back) return 1;

    // since xb is updated, we update d here.
    int i, c=0, k=0;
    int n = param[0], p = param[1], ns = param[3], set[ns];
    double value, *x_i, *x_j, ell[ns], lossafter;

    // calculate the first order Taylor's Formula
    der(x, y, d, meanx, sdx, model, param, residual, xb, isq, stautus, active,
        hierarchy);

    // saves active set's first order Taylor's Formula in ell and sort it.
    for (i = 0; i < ns; ++i) {
        c           = active[i];
        set[i]      = i;
        betanew[i]  = beta[i];
        ell[i]      = d[c] / weight[c];
        if (beta[i] > 0) ell[i] *= -1.0;
    }
    HeapSortDouble(ell, set, ns);

    // try a backward step by sort.
    ns--;
    while(ns >= 0)
    {
        c = set[ns];
        k = active[c];
        if (k >= p) break;
        if (fabs(beta[c]) > eps/weight[k] + xi) break;
        if (child[k] == 0) break;
        if (!hierarchy) {
            if (b_qlf(k, param, active, parent, child)) break;
        }
        ns--;
    }
    if (ns == -1) return 1; // none of active element can backward
    value = eps/weight[k];
    if (beta[c] > 0) value *= -1.0;
    betanew[c] += value;

    // calculate the object function and determine if we take a backward step
    value /= sdx[k];
    if (k < p)
    {
        x_i = x + k*n;
        for (i = 0; i < n; ++i)
            residual[i] = xb[i] + (x_i[i] - meanx[k])*value;
    } else {
        x_i = x + parent[2*c]*n;
        x_j = x + parent[2*c+1]*n;
        for (i = 0; i < n; ++i)
            residual[i] = xb[i] + (x_i[i]*x_j[i] - meanx[k]) * value;
    }
    lossafter = Loss(y, residual, model, param, stautus);

    if (lossafter - *losses - (*lambda)*eps < -1.0*xi)
    {
        // adopt a backward step;
        // update xb, lambda, beta, losses, direction, bic, df, child
        // param[3] (nonzero elements for step t), param[4] (nonzero main effect)
        losses[1]    = lossafter;
        direction[1] = -1;
        df[1]        = df[0];
        lambda[1]    = lambda[0];
        for (i = 0; i < n; ++i)
            xb[i] = residual[i];
        // test if vanish
        if (fabs(betanew[c]) < xi) {
            for (i = c; i < param[3]-1; ++i)
            {
                active[i]     = active[i+1];
                betanew[i]    = betanew[i+1];
                parent[2*i]   = parent[2*i+2];
                parent[2*i+1] = parent[2*i+3];
            }
            param[3]--;
            if (k < p) {
                param[4]--;
            } else {
                child[ parent[2*c] ]--;
                child[ parent[2*c+1] ]--;
            }
            df[1]--;
        }
        bic[1]       = calculate_bic(losses+1, param, model, gamma);
        return 0;
    }

    betanew[c] = beta[c];
    return 1;
}

// used in forward, to find the forward coordinate
void compare(double *temp, double value, int *k, int h12, int *k2, int t, 
    int *parenti, int i, int *parentj, int j)
{
    if (value > *temp)
    {
        *temp    = value;
        *k       = h12;
        *k2      = t;
        *parenti = i;
        *parentj = j;
    }
}

// take an forward step
void Forward(x, y, xb, beta, d, weight, param, model, eps, lambda, losses, 
    parent, xi, back, hierarchy, direction, active, bic, df, child, meanx, 
    sdx, residual, isq, stautus, gamma)
double *x, *y, *xb, *beta, *d, *weight, eps, *lambda, *losses, xi, *bic;
double *meanx, *sdx, *residual, gamma;
int *param, model, *parent, back, hierarchy, *direction, *active;
int *df, *child, isq, *stautus;
{
    // d doesn't need initial value
    int i, j, h12, k = -1, parenti = -1, parentj = -1, k2 = -1;
    int n = param[0], p = param[1], ns = param[3], nm = param[4];
    double value, value2, value3, temp, *x_i, *x_j;
    int r, ell = 0;

    // d has calculated in Backward.
    // beta has been assigned in backward
    if (!back) {
        der(x, y, d, meanx, sdx, model, param, residual, xb, isq, stautus, 
            active, hierarchy);
        for (i = 0; i < ns; ++i) beta[i] = beta[i-ns];
    }

    // find a forward direction
    temp = 0.0;
    for (i = 0; i < p; ++i)
    {  // main effects term
        value = fabs(d[i])/weight[i];
        if (value > temp)
        {
            temp = value;
            k = i;
        }
    }
    if (hierarchy)
    {
        for (i = 0; i <= active[nm-1]; ++i){
            x_i = x + i*n;
            if (i < active[ell])
            {
                r = (2*(p+isq-1)-i)*(i+1)/2;
                for (j = ell; j < nm; ++j) {
                    x_j = x + active[j]*n;
                    h12 = r + active[j];
                    value  = fabs(d[h12])/weight[h12]/2;
                    value2 = fabs(d[i])/weight[i]/2;
                    compare(&temp, value+value2, &k, h12, &k2, i, &parenti, i, &parentj, active[j]);
                }
            } else {
                j = i-isq+1;
                h12 = (2*(p+isq)-i)*(i+1)/2-isq;
                r = ell-isq+1;
                for (; j < p; ++j) {
                    x_j = x + j*n;
                    if ( (r < nm) && (j == active[r]) )
                    {
                        value = fabs(d[h12])/weight[h12];
                        compare(&temp, value, &k, h12, &k2, -1, &parenti, i, &parentj, j);
                        r++;
                    } else {
                        value = fabs(d[h12])/weight[h12]/2;
                        value2 = fabs(d[j])/weight[j]/2;
                        compare(&temp, value+value2, &k, h12, &k2, j, &parenti, i, &parentj, j);
                    }
                    h12++;
                }
                ell++;
            }
        }
    } else {
        // Weak hierarchy
        for (i = 0, x_i = x; i < p; ++i, x_i += n){
            //x_i  = x + i*n;
            j = i-isq+1;
            x_j = x + j*n;
            h12 = (2*(p+isq)-i)*(i+1)/2-isq;
            r = ell-isq+1;
            if ( (ell < nm) && (i == active[ell]) ) {
                ell++;
                for (; j < p; ++j, x_j += n) {
                    value = fabs(d[h12])/weight[h12];
                    compare(&temp, value, &k, h12, &k2, -1, &parenti, i, &parentj, j);
                    if ( (r < nm) && (j == active[r]) ) {
                        r++;
                    } else {
                        value2 = fabs(d[j])/weight[j];
                        compare(&temp, (value+value2)/2, &k, h12, &k2, j, &parenti, i, &parentj, j);
                    }
                    h12++;
                }
            } else {
                value3 = fabs(d[i])/weight[i];
                for (; j < p; ++j, x_j += n) {
                    value = fabs(d[h12])/weight[h12];
                    if ( (r < nm) && (j == active[r]) ) {
                        compare(&temp, value, &k, h12, &k2, -1, &parenti, i, &parentj, j);
                        r++;
                    } else {
                        value2 = fabs(d[j])/weight[j];
                        compare(&temp, (value+value2)/2, &k, h12, &k2, j, &parenti, i, &parentj, j);
                        compare(&temp, (value+value3)/2, &k, h12, &k2, i, &parenti, i, &parentj, j);
                    }
                    h12++;
                }
            }
        }
    }

    // calculate increment, update beta, active, df, param[3], param[4], child
    if (k2 != -1)
    {
        value = eps/weight[k];
        if (d[k] > 0.0) value *= -1.0;
        // now, k >= p, k2 < p, k,k2 not in active
        value2 = eps/weight[k2];
        if (d[k2] > 0.0) value2 *= -1.0;

        df[1]=df[0]+2;
        param[3]+=2;
        param[4]++;


        for (i = nm; i < ns; ++i)
        {
            if (active[i] < k) continue;
            break;
        }
        for (j = ns+1; j > i+1; --j)
        {
            active[j]     = active[j-2];
            beta[j]       = beta[j-2];
            parent[2*j]   = parent[2*j-4];
            parent[2*j+1] = parent[2*j-3];
        }
        active[i+1]     = k;
        beta[i+1]       = value;
        parent[2*i+2]   = parenti;
        parent[2*i+3]   = parentj;
        child[parenti]++;
        child[parentj]++;
        for (j = i; j > nm; --j)
        {
            active[j]     = active[j-1];
            beta[j]       = beta[j-1];
            parent[2*j]   = parent[2*j-2];
            parent[2*j+1] = parent[2*j-1];
        }

        for (i = 0; i < nm; ++i)
        {
            if (active[i] < k2) continue;
            break; 
        }
        for (j = nm; j > i; --j)
        {
            active[j]     = active[j-1];
            beta[j]       = beta[j-1];
        }
        active[i] = k2;
        beta[i]   = value2;

        // update xb
        value  /= sdx[k];
        value2 /= sdx[k2];
        x_i = x + k2*n;
        for (i = 0; i < n; ++i)
            xb[i] += (x_i[i] - meanx[k2])*value2;
        x_i = x + parenti*n;
        x_j = x + parentj*n;
        for (i = 0; i < n; ++i)
            xb[i] += (x_i[i]*x_j[i] - meanx[k])*value;
        direction[1] = 2;
    } else {
        value = eps/weight[k];
        if (d[k] > 0.0) value *= -1.0;

        if (k > active[ns-1])
        {
            active[ns] = k;
            beta[ns] = value;
            df[1] = df[0]+1;
            param[3]++;
            if (k < p)
            {
                param[4]++;
            } else {
                child[parenti]++;
                child[parentj]++;
                parent[2*ns] = parenti;
                parent[2*ns+1] = parentj;
            }
        } else {
            for (i = 0; i < ns; ++i)
            {
                if (active[i] < k) continue;
                if (active[i] == k)
                {
                    beta[i] += value;
                    df[1] = df[0];
                } else {
                    for (j = ns; j > i; --j)
                    {
                        active[j] = active[j-1];
                        beta[j] = beta[j-1];
                        parent[2*j] = parent[2*j-2];
                        parent[2*j+1] = parent[2*j-1];
                    }
                    active[i] = k;
                    beta[i] = value;
                    param[3]++;
                    df[1] = df[0]+1;
                    if (k < p)
                    {
                        param[4]++;
                    } else {
                        child[parenti]++;
                        child[parentj]++;
                        parent[2*i] = parenti;
                        parent[2*i+1] = parentj;
                    }
                }
                break;
            }
        }

        // update xb
        value /= sdx[k];
        if (k < p)
        {
            x_i = x + k*n;
            for (i = 0; i < n; ++i)
                xb[i] += (x_i[i] - meanx[k])*value;
        } else {
            x_i = x + parenti*n;
            x_j = x + parentj*n;
            for (i = 0; i < n; ++i)
                xb[i] += (x_i[i]*x_j[i] - meanx[k])*value;
        }
        direction[1] = 1;
    }

    // update lambda, losses, direction, bic, df
    losses[1]    = Loss(y, xb, model, param, stautus);
    //direction[1] = 1;
    bic[1]       = calculate_bic(losses+1, param, model, gamma);
    if (k2 != -1) {
        temp         = (losses[0] - losses[1])/(eps*2.0);
    } else {
        temp         = (losses[0] - losses[1])/eps;
    }
    lambda[1]    = temp < *lambda ? temp : *lambda;
}

int HFabs(x, y, xb, beta, d, weight, param, model, eps, lambda, losses, parent, 
    xi, back, hierarchy, direction, active, bic, df, child, meanx, sdx, residual, 
    isq, stautus, sparse_i, sparse_j, iter, stoping, lam_m, max_s, gamma)
double *x, *y, *xb, *beta, *d, *weight, eps, *lambda, *losses, xi, *bic;
double *meanx, *sdx, *residual, lam_m, gamma;
int *param, model, *parent, back, hierarchy, *direction, *active, stoping;
int *df, *child, isq, *stautus, *sparse_i, *sparse_j, *iter, max_s;
{
    int i, k;
    //setup_parent(parent, param, isq);
    Pre_work(x, y, weight, param, meanx, sdx, isq);

    // step 1: initial step (forward)
    Initial(x, y, xb, beta, d, weight, param, model, eps, lambda, losses, direction, 
        active, bic, df, meanx, sdx, residual, isq, stautus, hierarchy, gamma);
    int tt_act_b = 0;
    int tt_act_a = 0;
    usi(sparse_i, sparse_j, &tt_act_b, &tt_act_a, active, param[3], 0);

    // step 2: forward and backward
    for (i = 0; i < *iter-1; ++i)
    {
        k = Backward(x, y, xb, beta+tt_act_b, d, beta+tt_act_a, weight, param, 
            model, eps, lambda+i, losses+i, parent, xi, back, hierarchy, active, 
            direction+i, bic+i, df+i, child, meanx, sdx, residual, isq, stautus, gamma);
        if (k)
        {
            Forward(x, y, xb, beta+tt_act_a, d, weight, param, model, eps, 
                lambda+i, losses+i, parent, xi, back, hierarchy, direction+i, 
                active, bic+i, df+i, child, meanx, sdx, residual, isq, stautus, gamma);
        }
        usi(sparse_i, sparse_j, &tt_act_b, &tt_act_a, active, param[3], i+1);

        if ( stoping && (lambda[i+1] <= lambda[0] * lam_m) ) {
            *iter = i+2;
            if (lambda[i+1] < 0) 
            {
                (*iter)--;
                tt_act_a -= param[3];
            }
            break;
        }
        if (param[3] > max_s) {
            //Rprintf("Warning! Max nonzero number is larger than predetermined threshold. Program ended early.\n");
            *iter = i+2;
            break;
        }
        if (i == *iter-2) {
            Rprintf("Solution path unfinished, more iterations are needed.\n");
            break;
        }
    }
    return tt_act_a;
}

void Pre_work_GE(double *x, double *z, double *y, double *weight, int *param, 
    double *meanx, double *sdx)
{
    // y has been centralized, while x hasn't.
    int n = param[0];
    int px = param[1], pz = param[5];
    int i, j, k, c;
    double *x_i, *x_j, temp, var, cor;

    // mean and sd of X
    for (c = 0; c < px; ++c)
    {
        cor = 0.0;
        var = 0.0;
        x_i = x + c*n;
        for (k = 0; k < n; ++k)
        {
            cor += x_i[k];
            var += x_i[k]*x_i[k];
        }
        meanx[c] = cor/n;
        sdx[c] = sqrt( (var - n*meanx[c]*meanx[c])/(n-1) );
    }

    // mean and sd of Z
    for (; c < px+pz; ++c)
    {
        cor = 0.0;
        var = 0.0;
        x_i = z + (c-px)*n;
        for (k = 0; k < n; ++k)
        {
            cor += x_i[k];
            var += x_i[k]*x_i[k];
        }
        meanx[c] = cor/n;
        sdx[c] = sqrt( (var - n*meanx[c]*meanx[c])/(n-1) );
    }

    for (i = 0; i < px; ++i)
    {
        x_i = x + i*n;
        for (j = 0; j < pz; ++j)
        {
            x_j = z + j*n;
            cor = 0.0;
            var = 0.0;
            for (k = 0; k < n; ++k) 
            {
                temp = x_i[k]*x_j[k];
                cor += temp;
                var += temp*temp;
            }
            meanx[c] = cor/n;
            sdx[c] = sqrt( (var - n*meanx[c]*meanx[c])/(n-1) );
            c++;
        }
    }
}

// backward qualify for weak hierarchy, judge if a main effect k can vanish 
int b_qlf_GE(int k, int *param, int *active, int *parent, int *child)
{
    int i, j, prt, loc;
    int px = param[1], ns = param[3];
    int dfx = param[4], dfz = param[6], c = child[k];

    if (k < px) {
        loc = 0;
    } else {
        loc = 1;
    }

    for (i = dfx+dfz; i < ns; ++i)
    {
        prt = parent[2*i+loc];
        if ( prt == k ) {
            c--;
            prt = parent[2*i+1-loc];
            for (j = 0; j < dfx+dfz; ++j)
            {
                if (active[j] < prt) continue;
                if (active[j] > prt) return 0;
                break;
            }
        }
        if (c == 0) break;
    }
    return 1;
}

void der_GE(double *x, double *z, double *y, double *d, double *meanx, double *sdx, int model, 
    int *param, double *residual, double *xb, int *stautus, int hierarchy, int *active) 
{
    // calculate the 1st order Taylor Formula.
    int i, j, k, c, ellx=0, ellz=0;
    // Param = c(n, px, q, df, dfx, pz, dfz)
    int n = param[0], px = param[1], pz = param[5];
    int dfx = param[4], dfz = param[6];
    double *x_i, *x_j, temp;

    if(model == 1){
        temp = 0.0;
        for (i = 0; i < n; ++i) {
            residual[i] = y[i] - xb[i];
            temp += residual[i];
        }
        for (c = 0; c < px; ++c) {
            x_i = x + c*n;
            d[c] = 0.0;
            for (j = 0; j < n; ++j)
                d[c] -= x_i[j] * residual[j];
            d[c] = (d[c] + meanx[c]*temp)/(n*sdx[c]);
        }
        for (; c < px+pz; ++c) {
            x_i = z + (c-px)*n;
            d[c] = 0.0;
            for (j = 0; j < n; ++j)
                d[c] -= x_i[j] * residual[j];
            d[c] = (d[c] + meanx[c]*temp)/(n*sdx[c]);
        }
        if (hierarchy) {
            // Strong hierarchy
            for (i = 0; i < px; ++i)
            {
                x_i = x + i*n;
                if ( (ellx == dfx) || (i < active[ellx]) )
                {
                    //ellz = px + pz + i*pz - px;
                    ellz = (i+1)*pz;
                    for (j = dfx; j < dfz+dfx; ++j)
                    {
                        x_j = z + (active[j]-px)*n;
                        c = ellz + active[j];
                        d[c] = 0.0;
                        for (k = 0; k < n; ++k)
                            d[c] -= x_i[k]*x_j[k] * residual[k];
                        d[c] = (d[c] + meanx[c]*temp)/(n*sdx[c]);
                    }
                } else {
                    c = px + pz + i*pz;
                    for (j = 0; j < pz; ++j)
                    {
                        x_j = z + j*n;
                        d[c] = 0.0;
                        for (k = 0; k < n; ++k)
                            d[c] -= x_i[k]*x_j[k] * residual[k];
                        d[c] = (d[c] + meanx[c]*temp)/(n*sdx[c]);
                        c++;
                    }
                    ellx++;
                }
            }
        } else {
            // Weak hierarchy
            for (i = 0; i < px; ++i)
            {
                x_i = x + i*n;
                for (j = 0; j < pz; ++j)
                {
                    x_j = z + j*n;
                    d[c] = 0.0;
                    for (k = 0; k < n; ++k)
                        d[c] -= x_i[k]*x_j[k] * residual[k];
                    d[c] = (d[c] + meanx[c]*temp)/(n*sdx[c]);
                    c++;
                }
            }
            
        }
    } else if(model == 2) {
        double theta = 0.0;
        int ix, iz;

        for (i = 0; i < param[2]; ++i) d[i] = 0.0;

        for (i = n-1; i >= 0; --i) {
            theta += exp(xb[i]);
            if (stautus[i]) {
                for (j = i; j < n; ++j) {
                    x_i = x;
                    temp = exp(xb[j])/theta;
                    for (c = 0; c < px; ++c) {
                        d[c] -= temp * (x_i[i]-x_i[j]) / sdx[c];
                        x_i += n;
                    }
                    x_i = z;
                    for (; c < px+pz; ++c) {
                        d[c] -= temp * (x_i[i]-x_i[j]) / sdx[c];
                        x_i += n;
                    }

                    if (hierarchy)
                    {
                        // Strong hierarchy
                        ellx = 0;
                        for (ix = 0; ix < px; ++ix)
                        {
                            x_i = x + ix*n;
                            if ( (ellx == dfx) || (ix < active[ellx]) )
                            {
                                ellz = pz + ix*pz;
                                for (iz = dfx; iz < dfz+dfx; ++iz)
                                {
                                    x_j   = z + (active[iz]-px)*n;
                                    c     = ellz + active[iz];
                                    d[c] -= temp*(x_i[i]*x_j[i]-x_i[j]*x_j[j]) / sdx[c];
                                }
                            } else {
                                c   = px + pz + ix*pz;
                                x_j = z;
                                for (iz = 0; iz < pz; ++iz)
                                {
                                    d[c] -= temp*(x_i[i]*x_j[i]-x_i[j]*x_j[j]) / sdx[c];
                                    c++;
                                    x_j += n;
                                }
                                ellx++;
                            }
                        }
                    } else {
                        x_i = x;
                        for (k = 0; k < px; ++k){
                            x_j = z;
                            for (iz = 0; iz < pz; ++iz) {
                                d[c] -= temp*(x_i[i]*x_j[i]-x_i[j]*x_j[j]) / sdx[c];
                                c++;
                                x_j += n;
                            }
                            x_i += n;
                        }
                    }
                }
            }
        }
    }
}

// take an initial step
void Initial_GE(x, z, y, xb, beta, d, w, param, model, eps, lambda, losses, direction, 
    active, bic, df, meanx, sdx, residual, stautus, hierarchy, gamma)
double *x, *z, *y, *xb, *beta, *d, *w, eps, *lambda, *losses, *bic;
double *meanx, *sdx, *residual, gamma;
int *param, model, *direction, *active, *df, *stautus, hierarchy;
{
    // d doesn't need initial value
    int i, k=0;
    int n = param[0], px = param[1], pz = param[5];
    double value, temp, *xi;

    // calculate the derivative
    der_GE(x, z, y, d, meanx, sdx, model, param, residual, xb, stautus, hierarchy, active);

    // find a forward direction
    temp = 0.0;
    for (i = 0; i < px+pz; ++i){
        value = fabs(d[i])/w[i];
        if (value > temp){
            temp = value;
            k = i;
        }
    }

    // calculate increment, update xb, lambda, beta, losses, direction, bic, df
    // param[3] (nonzero elements for step t), param[4] (nonzero main effect)
    value = eps/w[k];
    if (d[k] > 0.0) value *= -1.0;
    *beta = value;

    temp = Loss(y, xb, model, param, stautus);
    if (k < px) {
        xi = x + k*n;
        param[4] = 1;
    } else {
        xi = z + (k-px)*n;
        param[6] = 1;
    }
    value /= sdx[k];
    for (i = 0; i < n; ++i)
        xb[i] += (xi[i]-meanx[k])*value;
    *losses = Loss(y, xb, model, param, stautus);
    *lambda = (temp - *losses)/eps;

    *direction = 1;
    param[3]   = 1;
    *active    = k;
    *bic       = calculate_bic(losses, param, model, gamma);
    *df        = 1;
}

// try a backward step
int Backward_GE(x, z, y, xb, beta, d, betanew, weight, param, model, eps, lambda, 
    losses, parent, xi, back, hierarchy, active, direction, bic, df, 
    child, meanx, sdx, residual, stautus, gamma)
double *x, *z, *y, *xb, *beta, *d, *betanew, *weight, eps, *lambda, *losses, xi;
double *bic, *meanx, *sdx, *residual, gamma;
int *param, model, *parent, back, hierarchy, *active, *direction;
int *df, *child, *stautus;
{
    if (!back) return 1;

    // since xb is updated, we update d here.
    int i, c=0, k=0;
    int n = param[0], px = param[1], pz = param[5], p = px+pz;
    int ns = param[3], set[ns];
    double value, *x_i, *x_j, ell[ns], lossafter;

    // calculate the first order Taylor's Formula
    der_GE(x, z, y, d, meanx, sdx, model, param, residual, xb, stautus, hierarchy, active);

    // saves active set's first order Taylor's Formula in ell and sort it.
    for (i = 0; i < ns; ++i) {
        c           = active[i];
        set[i]      = i;
        betanew[i]  = beta[i];
        ell[i]      = d[c] / weight[c];
        if (beta[i] > 0) ell[i] *= -1.0;
    }
    HeapSortDouble(ell, set, ns);

    // try a backward step by sort.
    ns--;
    while(ns >= 0)
    {
        c = set[ns];
        k = active[c];
        if (k >= p) break;
        if (fabs(beta[c]) > eps/weight[k] + xi) break;
        if (child[k] == 0) break;
        if (!hierarchy) {
            if (b_qlf_GE(k, param, active, parent, child)) break;
        }
        ns--;
    }
    //if (ns == -1) return 1; // none of active element can 
    value = eps/weight[k];
    if (beta[c] > 0) value *= -1.0;
    betanew[c] += value;

    // calculate the object function and determine if we take a backward step
    value /= sdx[k];
    if (k < p)
    {
        if (k < px) {
            x_i = x + k*n;
        } else {
            x_i = z + (k - px)*n;
        }
        for (i = 0; i < n; ++i)
            residual[i] = xb[i] + (x_i[i] - meanx[k])*value;
    } else {
        x_i = x + parent[2*c]*n;
        x_j = z + (parent[2*c+1]-px)*n;
        for (i = 0; i < n; ++i)
            residual[i] = xb[i] + (x_i[i]*x_j[i] - meanx[k]) * value;
    }
    lossafter = Loss(y, residual, model, param, stautus);

    if (lossafter - *losses - (*lambda)*eps < -1.0*xi)
    {
        // adopt a backward step;
        // update xb, lambda, beta, losses, direction, bic, df, child
        // param[3] (nonzero elements for step t), param[4] (nonzero main effect)
        losses[1]    = lossafter;
        direction[1] = -1;
        df[1]        = df[0];
        lambda[1]    = lambda[0];
        for (i = 0; i < n; ++i)
            xb[i] = residual[i];
        // test if vanish
        if (fabs(betanew[c]) < xi) {
            for (i = c; i < param[3]-1; ++i)
            {
                active[i]     = active[i+1];
                betanew[i]    = betanew[i+1];
                parent[2*i]   = parent[2*i+2];
                parent[2*i+1] = parent[2*i+3];
            }
            param[3]--;
            if (k < p) {
                if (k < px) {
                    param[4]--;
                } else {
                    param[6]--;
                }
            } else {
                child[ parent[2*c] ]--;
                child[ parent[2*c+1] ]--;
            }
            df[1]--;
        }
        bic[1]       = calculate_bic(losses+1, param, model, gamma);
        return 0;
    }

    betanew[c] = beta[c];
    return 1;
}

// take an forward step
void Forward_GE(x, z, y, xb, beta, d, weight, param, model, eps, lambda, losses, 
    parent, xi, back, hierarchy, direction, active, bic, df, child, meanx, 
    sdx, residual, stautus, gamma)
double *x, *z, *y, *xb, *beta, *d, *weight, eps, *lambda, *losses, xi, *bic;
double *meanx, *sdx, *residual, gamma;
int *param, model, *parent, back, hierarchy, *direction, *active;
int *df, *child, *stautus;
{
    // d doesn't need initial value
    int i, j, h12, ellz = 0, ellx = 0, k = -1, k2 = -1;
    int n = param[0], px = param[1], pz = param[5], p = px+pz;
    int ns = param[3], dfx = param[4], dfz = param[6];
    int parenti = -1, parentj = -1;
    double value, value2, temp, *x_i, *x_j;

    // d has calculated in Backward.
    // beta has been assigned in backward
    if (!back) {
        der_GE(x, z, y, d, meanx, sdx, model, param, residual, xb, stautus, hierarchy, active);
        for (i = 0; i < ns; ++i) beta[i] = beta[i-ns];
    }

    // find a forward direction
    temp = 0.0;
    for (i = 0; i < p; ++i)
    {   // main effects term
        value = fabs(d[i])/weight[i];
        if (value > temp)
        {
            temp = value;
            k = i;
        }
    }
    
    if (hierarchy)
    {
        for (i = 0; i < px; ++i){
            x_i = x + i*n;
            if ( (ellx == dfx) || (i < active[ellx]) )
            {
                // x is not an active coordinate
                ellz = pz + i*pz;
                value2 = fabs(d[i])/weight[i]/2;
                for (j = dfx; j < dfz+dfx; ++j) {
                    // z is an active coordinate
                    x_j = z + (active[j]-px)*n;
                    h12 = ellz + active[j];
                    value  = fabs(d[h12])/weight[h12]/2;
                    compare(&temp, value+value2, &k, h12, &k2, i, &parenti, i, &parentj, active[j]);
                }
            } else {
                h12 = p + i*pz;
                ellz = dfx;
                for (j = px; j < p; ++j) {
                    x_j = z + (j-px)*n;
                    if ( (ellz == dfz+dfx) || (j < active[ellz]) )
                    {
                        // z is not an active coordinate
                        value  = fabs(d[h12])/weight[h12]/2;
                        value2 = fabs(d[j])/weight[j]/2;
                        compare(&temp, value+value2, &k, h12, &k2, j, &parenti, i, &parentj, j);
                    } else {
                        value = fabs(d[h12])/weight[h12];
                        compare(&temp, value, &k, h12, &k2, -1, &parenti, i, &parentj, j);
                        ellz++;
                    }
                    h12++;
                }
                ellx++;
            }
        }
    } else {
        // Weak hierarchy
        for (i = 0; i < px; ++i){
            x_i  = x + i*n;
            h12  = p + i*pz;
            ellz = dfx;
            if ( (ellx < dfx) && (i == active[ellx]) ) {
                // x is an active coordinate
                ellx++;
                for (j = px; j < p; ++j) {
                    x_j = z + (j-px)*n;
                    value = fabs(d[h12])/weight[h12];
                    compare(&temp, value, &k, h12, &k2, -1, &parenti, i, &parentj, j);
                    if ( (ellz == dfx+dfz) || (j < active[ellz]) )
                    {
                        // z is not an active coordinate
                        value2 = fabs(d[j])/weight[j];
                        compare(&temp, (value+value2)/2, &k, h12, &k2, j, &parenti, i, &parentj, j);
                    } else {
                        ellz++;
                    }
                    h12++;
                }
            } else {
                // x is not an active coordinate
                for (j = px; j < p; ++j) {
                    x_j = z + (j-px)*n;
                    value = fabs(d[h12])/weight[h12];
                    if ( (ellz == dfx+dfz) || (j < active[ellz]) )
                    {
                        // z is not an active coordinate
                        value2 = fabs(d[j])/weight[j];
                        compare(&temp, (value+value2)/2, &k, h12, &k2, j, &parenti, i, &parentj, j);
                        value2 = fabs(d[i])/weight[i];
                        compare(&temp, (value+value2)/2, &k, h12, &k2, i, &parenti, i, &parentj, j);
                    } else {
                        compare(&temp, value, &k, h12, &k2, -1, &parenti, i, &parentj, j);
                        ellz++;
                    }
                    h12++;
                }
            }
        }
    }

    // calculate increment, update beta, active, df, param[3], param[4], child
    if (k2 != -1)
    {
        value = eps/weight[k];
        if (d[k] > 0.0) value *= -1.0;
        // now, k >= p, k2 < p, k,k2 not in active
        value2 = eps/weight[k2];
        if (d[k2] > 0.0) value2 *= -1.0;

        df[1]=df[0]+2;
        param[3]+=2;

        for (i = dfx+dfz; i < ns; ++i)
        {
            if (active[i] < k) continue;
            break;
        }
        // active[i] > k
        for (j = ns+1; j > i+1; --j)
        {
            active[j]     = active[j-2];
            beta[j]       = beta[j-2];
            parent[2*j]   = parent[2*j-4];
            parent[2*j+1] = parent[2*j-3];
        }
        active[i+1]     = k;
        beta[i+1]       = value;
        parent[2*i+2]   = parenti;
        parent[2*i+3]   = parentj;
        child[parenti]++;
        child[parentj]++;

        for (j = 0; j < dfx+dfz; ++j)
        {
            if (active[j] < k2) continue;
            break; 
        }
        // active[j] > k2
        for (h12 = i; h12 > j; --h12)
        {
            active[h12]     = active[h12-1];
            beta[h12]       = beta[h12-1];
            parent[2*h12]   = parent[2*h12-2];
            parent[2*h12+1] = parent[2*h12-1];
        }
        active[j] = k2;
        beta[j]   = value2;
        
        if (k2 < px) {
            param[4]++;
            x_i = x + k2*n;
        } else {
            param[6]++;
            x_i = z + (k2-px)*n;
        }

        // update xb
        value2 /= sdx[k2];
        for (i = 0; i < n; ++i)
            xb[i] += (x_i[i] - meanx[k2])*value2;
        value  /= sdx[k];
        x_i = x + parenti*n;
        x_j = z + (parentj-px)*n;
        for (i = 0; i < n; ++i)
            xb[i] += (x_i[i]*x_j[i] - meanx[k])*value;
        direction[1] = 2;
    } else {
        value = eps/weight[k];
        if (d[k] > 0.0) value *= -1.0;

        for (i = 0; i < ns; ++i)
        {
            if (active[i] < k) continue;
            if (active[i] == k)
            {
                beta[i] += value;
                df[1] = df[0];
                ellz = -1;
            }
            break;
        }

        if (ellz != -1)
        {
            for (j = ns; j > i; --j)
            {
                active[j]     = active[j-1];
                beta[j]       = beta[j-1];
                parent[2*j]   = parent[2*j-2];
                parent[2*j+1] = parent[2*j-1];
            }
            active[i] = k;
            beta[i]   = value;
            df[1]     = df[0]+1;
            param[3]++;
            if (k < p)
            {
                if (k < px) {
                    param[4]++;
                } else {
                    param[6]++;
                }
            } else {
                child[parenti]++;
                child[parentj]++;
                parent[2*i]   = parenti;
                parent[2*i+1] = parentj;
            }
        }
        // update xb
        value /= sdx[k];
        if (k < p)
        {
            if (k < px) {
                x_i = x + k*n;
            } else {
                x_i = z + (k-px)*n;
            }
            for (i = 0; i < n; ++i)
                xb[i] += (x_i[i] - meanx[k])*value;
        } else {
            x_i = x + parenti*n;
            x_j = z + (parentj-px)*n;
            for (i = 0; i < n; ++i)
                xb[i] += (x_i[i]*x_j[i] - meanx[k])*value;
        }
        direction[1] = 1;
    }

    // update lambda, losses, direction, bic, df
    losses[1]    = Loss(y, xb, model, param, stautus);
    //direction[1] = 1;
    bic[1]       = calculate_bic(losses+1, param, model, gamma);
    if (k2 != -1) {
        temp         = (losses[0] - losses[1])/(eps*2.0);
    } else {
        temp         = (losses[0] - losses[1])/eps;
    }
    lambda[1]    = temp < *lambda ? temp : *lambda;
}

int HFabs_GE(x, z, y, xb, beta, d, weight, param, model, eps, lambda, losses, 
    parent, xi, back, hierarchy, direction, active, bic, df, child, meanx, 
    sdx, residual, stautus, sparse_i, sparse_j, iter, stoping, lam_m, max_s, gamma)
double *x, *z, *y, *xb, *beta, *d, *weight, eps, *lambda, *losses, xi, *bic;
double *meanx, *sdx, *residual, lam_m, gamma;
int *param, model, *parent, back, hierarchy, *direction, *active, *iter;
int *df, *child, *stautus, *sparse_i, *sparse_j, stoping, max_s;
{
    // Param = c(n, px, q, df, dfx, pz, dfz)
    int i, k;
    Pre_work_GE(x, z, y, weight, param, meanx, sdx);

    // step 1: initial step (forward)
    Initial_GE(x, z, y, xb, beta, d, weight, param, model, eps, lambda, losses, direction, 
        active, bic, df, meanx, sdx, residual, stautus, hierarchy, gamma);
    int tt_act_b = 0;
    int tt_act_a = 0;
    usi(sparse_i, sparse_j, &tt_act_b, &tt_act_a, active, param[3], 0);

    // step 2: forward and backward
    for (i = 0; i < *iter-1; ++i)
    {
        k = Backward_GE(x, z, y, xb, beta+tt_act_b, d, beta+tt_act_a, weight, param, 
            model, eps, lambda+i, losses+i, parent, xi, back, hierarchy, active, 
            direction+i, bic+i, df+i, child, meanx, sdx, residual, stautus, gamma);
        if (k)
        {
            Forward_GE(x, z, y, xb, beta+tt_act_a, d, weight, param, model, eps, 
                lambda+i, losses+i, parent, xi, back, hierarchy, direction+i, 
                active, bic+i, df+i, child, meanx, sdx, residual, stautus, gamma);
        }
        usi(sparse_i, sparse_j, &tt_act_b, &tt_act_a, active, param[3], i+1);

        if ( stoping && (lambda[i+1] <= lambda[0] * lam_m) ) {
            *iter = i+2;
            if (lambda[i+1] < 0) 
            {
                (*iter)--;
                tt_act_a -= param[3];
            }
            break;
        }
        if (param[3] > max_s) {
            //Rprintf("Warning! Max nonzero number is larger than predetermined threshold. Program ended early.\n");
            *iter = i+2;
            break;
        }
        if (i == *iter-2) {
            Rprintf("Solution path unfinished, more iterations are needed.\n");
            break;
        }
    }
    return tt_act_a;
}

SEXP Hierarchy_Fabs(SEXP X, SEXP Z, SEXP Y, SEXP Weight, SEXP Model, SEXP Epsilon, SEXP Lam_min, 
    SEXP Xi, SEXP Back, SEXP Stoping, SEXP Iter, SEXP Hierarchy, SEXP Param, 
    SEXP Max_S, SEXP MeanY, SEXP Isquadratic, SEXP Status, SEXP GE, SEXP Gamma)
{
    int i, n, p, pz, q;
    double *x, *z, *y, *weight, eps, lam_m, xi, *meanx, *sdx, gamma;
    int *param, model, back, stoping, iter, hierarchy, max_s, *stautus, isq, ge;
    const char *modeltype = CHAR(asChar(Model));
    const char *Constrain = CHAR(asChar(Hierarchy));

    if (strcmp(modeltype, "gaussian") == 0) model = 1;
    else if (strcmp(modeltype, "cox") == 0) model = 2;
    else model = 1;

    if (strcmp(Constrain, "strong") == 0) hierarchy = 1;
    else if (strcmp(Constrain, "weak") == 0) hierarchy = 0;
    else hierarchy = 1;

    param     = INTEGER(Param);
    n         = param[0];
    p         = param[1];
    pz        = param[5];
    q         = param[2];
    x         = REAL(X);     // X has been scaled in R
    z         = REAL(Z);
    y         = REAL(Y);
    weight    = REAL(Weight);
    eps       = REAL(Epsilon)[0];
    lam_m     = REAL(Lam_min)[0];
    xi        = REAL(Xi)[0];
    back      = INTEGER(Back)[0];
    stoping   = INTEGER(Stoping)[0];
    iter      = INTEGER(Iter)[0];
    max_s     = INTEGER(Max_S)[0];
    isq       = INTEGER(Isquadratic)[0];
    stautus   = INTEGER(Status);
    meanx     = (double*)malloc(sizeof(double)  *q);
    sdx       = (double*)malloc(sizeof(double)  *q);
    ge        = INTEGER(GE)[0];
    gamma     = REAL(Gamma)[0];

    double *d, *beta, *lambda, *xb, *bic, *losses, *residual;
    int *direction, *df, *parent, *child, *active;
    int *sparse_i, *sparse_j, tt_act_a;

    //beta      = (double*)calloc(iter*max_s, sizeof(double));
    beta      = (double*)malloc(sizeof(double)*iter*max_s);
    sparse_i  =    (int*)calloc(iter*max_s, sizeof(int));
    sparse_j  =    (int*)calloc(iter*max_s, sizeof(int));
    lambda    = (double*)malloc(sizeof(double)*iter);
    direction =    (int*)malloc(sizeof(int)   *iter);
    bic       = (double*)malloc(sizeof(double)*iter);
    losses    = (double*)malloc(sizeof(double)*iter);
    df        =    (int*)malloc(sizeof(int)   *iter);
    residual  = (double*)malloc(sizeof(double)  *n);
    d         = (double*)malloc(sizeof(double)  *q);// 1st order Taylor Formula.
    xb        = (double*)calloc(n, sizeof(double)); // x beta, needs initialized
    child     =    (int*)calloc(p+pz, sizeof(int));    // number of child.
    parent    =    (int*)malloc(sizeof(int)   *2*(max_s+2));
    active    =    (int*)calloc(max_s+2, sizeof(int));// needs initialized

    if (ge) {
        tt_act_a = HFabs_GE(x, z, y, xb, beta, d, weight, param, model, eps, 
            lambda, losses, parent, xi, back, hierarchy, direction, active, bic, 
            df, child, meanx, sdx, residual, stautus, sparse_i, 
            sparse_j, &iter, stoping, lam_m, max_s, gamma);
    } else {
        tt_act_a = HFabs(x, y, xb, beta, d, weight, param, model, eps, lambda, 
            losses, parent, xi, back, hierarchy, direction, active, bic, df, 
            child, meanx, sdx, residual, isq, stautus, sparse_i, 
            sparse_j, &iter, stoping, lam_m, max_s, gamma);
    }

    SEXP Beta, Lambda, Direction, Loops, BIC, Loss, Df, Indexi, Indexj;
    SEXP Intercept, Result, R_names;
    char *names[10] = {"beta", "lambda", "direction", "iter", "bic", "loss", 
    "df", "index_i", "index_j", "intercept"};
    PROTECT(Beta      = allocVector(REALSXP, tt_act_a));
    PROTECT(Indexi    = allocVector(INTSXP,  tt_act_a));
    PROTECT(Indexj    = allocVector(INTSXP,  tt_act_a));
    PROTECT(Lambda    = allocVector(REALSXP, iter));
    PROTECT(BIC       = allocVector(REALSXP, iter));
    PROTECT(Loss      = allocVector(REALSXP, iter));
    PROTECT(Intercept = allocVector(REALSXP, iter));
    PROTECT(Direction = allocVector(INTSXP,  iter));
    PROTECT(Df        = allocVector(INTSXP,  iter));
    PROTECT(Loops     = allocVector(INTSXP,  1));
    PROTECT(Result    = allocVector(VECSXP,  10));
    PROTECT(R_names   = allocVector(STRSXP,  10));

    unmormalized(beta, REAL(Beta), parent, sparse_i, sparse_j, iter, 
        REAL(Intercept), meanx, sdx, REAL(MeanY)[0], tt_act_a, model);

    for(i = 0; i < 10; ++i) SET_STRING_ELT(R_names, i,  mkChar(names[i]));
    INTEGER(Loops)[0] = iter;
    for (i = 0; i < tt_act_a; ++i) 
    {
        INTEGER(Indexi)[i] = sparse_i[i];
        INTEGER(Indexj)[i] = sparse_j[i];
    }
    for (i = 0; i < iter; ++i) 
    {
        REAL(BIC)[i]          = bic[i];
        REAL(Loss)[i]         = losses[i];
        REAL(Lambda)[i]       = lambda[i];
        INTEGER(Direction)[i] = direction[i];
        INTEGER(Df)[i]        = df[i];
    }
    
    SET_VECTOR_ELT(Result, 0, Beta);
    SET_VECTOR_ELT(Result, 1, Lambda);
    SET_VECTOR_ELT(Result, 2, Direction);
    SET_VECTOR_ELT(Result, 3, Loops); 
    SET_VECTOR_ELT(Result, 4, BIC);
    SET_VECTOR_ELT(Result, 5, Loss);  
    SET_VECTOR_ELT(Result, 6, Df);  
    SET_VECTOR_ELT(Result, 7, Indexi);  
    SET_VECTOR_ELT(Result, 8, Indexj); 
    SET_VECTOR_ELT(Result, 9, Intercept);    
    setAttrib(Result, R_NamesSymbol, R_names); 


    free(meanx);
    free(sdx);
    free(beta);
    free(sparse_i);
    free(sparse_j);
    free(lambda);
    free(direction);
    free(bic);
    free(losses);
    free(df);
    free(parent);
    free(residual);
    free(d);
    free(xb);
    free(child);
    free(active);

    UNPROTECT(12);
    return Result;
}

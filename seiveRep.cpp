//
// Created by advanced on 17-1-16.
//

#include "MOPSO.h"

void sortRepByAngle(myRep &rep, ANGLE *angle){
    int len = static_cast<int>(rep.size());
    myRep::iterator it ;
    POINT point[len+100];
    int k=0;
    for (it = rep.begin(); it != rep.end(); it++){
        point[k].x = it->Cost[0];
        point[k].y = it->Cost[1];
        point[k].id = k;
        k++;
    }
    pointQsort(0, len-1, point);
    for (int i=0; i<len; i++){
        angle[i].id = point[i].id;
    }
    for (int i=1; i<len-1; i++){
        getAngle(i, point, angle);
    }
    angle[0].value     = -PI;
    angle[len-1].value = -2*PI;
    angleQsort(0, len - 1, angle);

    for (int i = 0; i < len; i++){
        sortAns[i] = angle[i].id;
    }
}

void getAngle(int i, const POINT *point, ANGLE *angle){
    VECTOR AB(point[i-1].x-point[i].x, point[i-1].y-point[i].y);
    VECTOR AC(point[i+1].x-point[i].x, point[i+1].y-point[i].y);

    double DOT = vectorDot(AB, AC);
    double l1 = length(AB), l2 = length(AC);
    double cosTheta = DOT / (l1 * l2);
    angle[i].value = acos(cosTheta);
    if (pdPlusPi(point[i-1], point[i], point[i+1]))
        angle[i].value = 2*PI - angle[i].value;
}

bool pdPlusPi(POINT first, POINT second, POINT third) {
    VECTOR BA(first.x - second.x, first.y - second.y), BC(third.x - second.x, third.y - second.y);
    if (cross(BA,BC) < 0) return 1;
    else return 0;
}

double vectorDot(VECTOR u, VECTOR v){
    return u.x * v.x + u.y * v.y;
}

double cross(VECTOR u, VECTOR v){
    return u.x * v.y - v.x * u.y ;              //  x1*y2 - x2*y1
}

double length(VECTOR u){
    return sqrt(u.x*u.x + u.y*u.y);
}

void angleQsort(int l, int r, ANGLE *a){
    if (l>=r)
        return;

    int i=l, j=r;
    double holeValue=a[i].value;
    int holeId = a[i].id;
    while (i<j){
        while (i<j && a[j].value <= holeValue) j--;
        if (i<j){
            a[i] = a[j];
            i++;
        }
        while (i<j && a[i].value >= holeValue) i++;
        if (i<j){
            a[j]= a[i];
            j--;
        }
    }

    a[i].value = holeValue;
    a[i].id = holeId;
    angleQsort(l, i - 1, a);
    angleQsort(i + 1, r, a);
}

void pointQsort(int l, int r, POINT *a){
    if (l>=r)
        return;

    int i=l, j=r;
    POINT holeValue = a[i];
    while (i<j){
        while (i<j && a[j].x >= holeValue.x) j--;
        if (i<j){
            a[i] = a[j];
            i++;
        }
        while (i<j && a[i].x <= holeValue.x) i++;
        if (i<j){
            a[j] = a[i];
            j--;
        }
    }

    a[i] = holeValue;
    pointQsort(l  , i-1, a);
    pointQsort(i+1, r, a);
}

void sortRepByLambda(myRep &rep){
    LAMBDA lambda[nRep + 100];
    myRep::iterator it = rep.begin();
    int i = 0;
    while (it != rep.end()){
        double totA = 0;
        double lambda1, lambda2;
        lambda1 = rand() / double(RAND_MAX);
        lambda2 = rand() / double(RAND_MAX);

        for (int j = 0; j < lambdaLoopTimes; j++){
            double A = lambda1 * it->Cost[0] + lambda2 * it->Cost[1] + (1 - lambda1 - lambda2) * it->Cost[2];
            totA += A;
            lambda1 *= lambda1;
            lambda2 *= lambda2;
        }

        lambda[i].value = totA / lambdaLoopTimes;
        lambda[i].id = i;
        i++;
        it++;
    }

    int len = static_cast<int>(rep.size());
    lambdaQsort(0,  len- 1, lambda);
    for (int j = 0; j < len; j++){
        sortAns[j] = lambda[j].id;
    }
}

void lambdaQsort(int l, int r, LAMBDA *a){
    if (l>=r)
        return;

    int i=l, j=r;
    double holeValue=a[i].value;
    int holeId = a[i].id;
    while (i<j){
        while (i<j && a[j].value <= holeValue) j--;
        if (i<j){
            a[i] = a[j];
            i++;
        }
        while (i<j && a[i].value >= holeValue) i++;
        if (i<j){
            a[j]= a[i];
            j--;
        }
    }

    a[i].value = holeValue;
    a[i].id = holeId;
    angleQsort(l, i - 1, a);
    angleQsort(i + 1, r, a);
}

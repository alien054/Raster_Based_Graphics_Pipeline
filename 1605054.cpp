#include<bits/stdc++.h>
#include "bitmap_image.hpp"

using namespace std;

#define mat_size 4
#define pi (2*acos(0.0))

int screen_width,screen_height;
double left_lim_x,bottom_lim_y;
double right_lim_x,top_lim_y;
double front_z,rear_z;
double dx,dy;
double top_y,bottom_y;
double left_x,right_x;

class Point4D
{
    public:

    double x,y,z,w;

    Point4D(){}

    Point4D(double x,double y,double z,double w)
    {
        this->x = x;
        this->y = y;
        this->z = z;
        this->w = w;
    }

    Point4D normalize()
    {
        double norm = sqrt((this->x*this->x) + (this->y*this->y) + (this->z*this->z));
        if(norm==0.0) norm = 0.0000000000001;

        double x_norm = (this->x*1.0)/(norm*1.0);
        double y_norm = (this->y*1.0)/(norm*1.0);
        double z_norm = (this->z*1.0)/(norm*1.0);

        Point4D normalized(x_norm,y_norm,z_norm,1);
        return normalized;
    }

    Point4D scale()
    {
        return Point4D(this->x/this->w,this->y/this->w,this->z/this->w,this->w/this->w);
    }

    Point4D resultant(Point4D other)
    {
        double x,y,z;
        x = this->x + other.x;
        y = this->y + other.y;
        z = this->z + other.z;

        Point4D res(x,y,z,1);
        return res;
    }
    float dot(Point4D other)
    {
        return (this->x*other.x + this->y*other.y + this->z*other.z);
    }

    Point4D scalar_product(double m)
    {
        double x,y,z;
        x = m * this->x;
        y = m * this->y;
        z = m * this->z;

        Point4D product(x,y,z,1);
        return product;
    }

    Point4D cross(Point4D other)
    {
        double x,y,z;
        x = this->y * other.z - this->z * other.y;
        y = this->z * other.x - this->x * other.z;
        z = this->x * other.y - this->y * other.x;

        Point4D cross_vect(x,y,z,1);
        return cross_vect;
    }



    Point4D rodrigues(Point4D unit_vec,double angle)
    {
        angle = (pi*angle*1.0)/180.0;

        //cos(t) * x
        Point4D cos_product = unit_vec.scalar_product(cos(angle));

        //(1-cos(t)) * (a.x) * a
        double a_dot_x = this->dot(unit_vec);
        double m = (1-cos(angle)) * a_dot_x;
        Point4D middle_product = this->scalar_product(m);

        //sin(t) * (a X x)
        Point4D a_cross_x = this->cross(unit_vec);
        Point4D sin_product = a_cross_x.scalar_product(sin(angle));

        //cos(t) * x + (1-cos(t)) * (a.x) * a + sin(t) * (a X x)
        return cos_product.resultant(middle_product.resultant(sin_product));
    }

    bool isParallel(Point4D other)
    {
        return (abs(this->y-other.y)<(dy/2));
    }

    void print()
    {
        cout<<this->x<<" "<<this->y<<" "<<this->z<<endl;
    }
};

class Color
{
public:
    int r,g,b;

    Color()
    {
        this->r = rand()%256;
        this->g = rand()%256;
        this->b = rand()%256;
    }

    Color(int r,int g,int b)
    {
        this->r = r;
        this->g = g;
        this->b = b;
    }

    void print()
    {
        cout<<"R: "<<this->r<<"\tG: "<<this->g<<"\tB: "<<this->b<<endl;
    }
};

class Triangle
{
public:
    Point4D p1,p2,p3;
    Color c;

    Triangle() {}
    Triangle(Point4D p1,Point4D p2,Point4D p3)
    {
        this->p1 = p1;
        this->p2 = p2;
        this->p3 = p3;
    }

    void print()
    {
        this->p1.print();
        this->p2.print();
        this->p3.print();
        this->c.print();
        cout<<endl;
    }
};

class Matrix
{
public:
    double matrix[mat_size][mat_size];

    Matrix()
    {
        for(int i=0;i<mat_size;i++)
        {
            for(int j=0;j<mat_size;j++)
            {
                if(i==j) this->matrix[i][j] = 1.0;
                else this->matrix[i][j] = 0;
            }
        }
    }

    Matrix(vector<double> value)
    {
        if(value.size() != 16)
        {
            cout<<"4X4 matrix must have 16 values"<<endl;
            exit(0);
        }

        for(int i=0;i<mat_size;i++)
        {
            for(int j=0;j<mat_size;j++)
            {
                this->matrix[i][j] = value[i*4+j];
            }
        }

    }

    Matrix multiply(Matrix other)
    {
        vector<double> mat;
        double res[mat_size][mat_size];

        for (int i=0; i<mat_size; i++)
        {
            for (int j=0; j<mat_size; j++)
            {
                res[i][j] = 0;
                for (int k = 0; k < mat_size; k++)
                {
                    res[i][j] += this->matrix[i][k] * other.matrix[k][j];
                }
            }
        }

        for(int i=0;i<mat_size;i++)
        {
            for(int j=0;j<mat_size;j++)
            {
                mat.push_back(res[i][j]);
            }
        }
        Matrix product(mat);
        return product;
    }

    Point4D multiply_point(Point4D point)
    {
        double p[5];
        double result[5];

        p[0] = point.x;
        p[1] = point.y;
        p[2] = point.z;
        p[3] = point.w;

        for(int i=0;i<mat_size;i++) result[i] = 0;

        for(int i=0;i<mat_size;i++)
        {
            for(int j=0;j<mat_size;j++)
            {
                result[i] += this->matrix[i][j] * p[j];
            }
        }

        Point4D temp(result[0],result[1],result[2],result[3]);
        Point4D res = temp.scale();

        return res;
    }

    void print()
    {
        cout<<"========================"<<endl;
        for(int i=0;i<mat_size;i++)
        {
            for(int j=0;j<mat_size;j++)
            {
                cout<<this->matrix[i][j]<<" ";
            }
            cout<<endl;
        }
        cout<<"========================"<<endl;
    }
};


int main(int argc, const char** argv) {
    srand (time(0));
    freopen("./test/1/scene.txt", "rb", stdin);
    string line,temp;
    Point4D eye,look,up,perspective;
    vector<Triangle> triangles;
    double temp1,temp2,temp3,temp4;
    Matrix model;
    stack<Matrix> mat_stack;
    bool push_flag = false;


    mat_stack.push(model);

    cin>>temp1>>temp2>>temp3;
    eye = Point4D(temp1,temp2,temp3,1);
    //cout<<"eye: "<<endl;
    //eye.print();

    cin>>temp1>>temp2>>temp3;
    look = Point4D(temp1,temp2,temp3,1);
    //cout<<"look: "<<endl;
    //look.print();

    cin>>temp1>>temp2>>temp3;
    up = Point4D(temp1,temp2,temp3,1);
    //cout<<"up: "<<endl;
    //up.print();

    cin>>temp1>>temp2>>temp3>>temp4;
    perspective = Point4D(temp1,temp2,temp3,temp4);
    //cout<<"perspective: "<<endl;
    //perspective.print();

    while(getline(cin, line))
    {
        //cout<<line<<endl;
        //removing \r from input
        if(line[line.length()-1] == '\r')
            line.erase(line.size() - 1);

        if(line == "triangle")
        {
            double t1,t2,t3;
            Point4D points[4];
            Point4D newPoints[4];
            for(int i=0;i<3;i++)
            {
                cin>>t1>>t2>>t3;
                points[i] = Point4D(t1,t2,t3,1);
            }

            for(int i=0;i<3;i++)
            {
                newPoints[i] = mat_stack.top().multiply_point(points[i]);
            }

            Triangle new_triangle(newPoints[0],newPoints[1],newPoints[2]);
            //new_triangle.print();
            triangles.push_back(new_triangle);
        }

        else if(line == "translate")
        {
            double value[4];
            vector<double> mat;

            cin>>value[0]>>value[1]>>value[2];
            value[3] = 1;

            for(int i=0;i<mat_size;i++)
            {
                for(int j=0;j<mat_size;j++)
                {
                    if(i==j) mat.push_back(1);
                    else if(j==mat_size-1) mat.push_back(value[i]);
                    else mat.push_back(0);
                }
            }

            Matrix transMat(mat);
            //transMat.print();

            Matrix top = mat_stack.top();
            Matrix temp_mat = top.multiply(transMat);

            if(push_flag)
            {
                mat_stack.push(temp_mat);
                push_flag = false;
            }
            else
            {
                mat_stack.pop();
                mat_stack.push(temp_mat);
            }

        }

        else if(line == "scale")
        {
            double value[4];
            vector<double> mat;

            cin>>value[0]>>value[1]>>value[2];
            value[3] = 1;

            for(int i=0;i<mat_size;i++)
            {
                for(int j=0;j<mat_size;j++)
                {
                    if(i==j) mat.push_back(value[i]);
                    else mat.push_back(0);
                }
            }

            Matrix scaleMat(mat);
            //scaleMat.print();

            Matrix top = mat_stack.top();
            Matrix temp_mat = top.multiply(scaleMat);;

            if(push_flag)
            {
                mat_stack.push(temp_mat);
                push_flag = false;
            }
            else
            {
                mat_stack.pop();
                mat_stack.push(temp_mat);
            }
        }

        else if(line == "rotate")
        {
            double value[4];
            double angle;
            Point4D columns[4];
            vector<double> mat;

            cin>>angle>>value[0]>>value[1]>>value[2];

            Point4D a(value[0],value[1],value[2],1);

            Point4D i(1,0,0,1);
            Point4D j(0,1,0,1);
            Point4D k(0,0,1,1);

            Point4D a_norm = a.normalize();
            columns[0] = a_norm.rodrigues(i,angle);
            columns[1] = a_norm.rodrigues(j,angle);
            columns[2] = a_norm.rodrigues(k,angle);

            for(int i=0;i<mat_size;i++)
            {
                for(int j=0;j<mat_size;j++)
                {
                    if(i!=3 && j==3) mat.push_back(0);
                    else if(i==3 && j==3) mat.push_back(1);
                    else
                    {
                        if(i==0) mat.push_back(columns[j].x);
                        else if(i==1) mat.push_back(columns[j].y);
                        else if(i==2) mat.push_back(columns[j].z);
                        else mat.push_back(0);
                    }
                }
            }

            Matrix rotMat(mat);
            //rotMat.print();

            Matrix top = mat_stack.top();
            Matrix temp_mat = top.multiply(rotMat);

            if(push_flag)
            {
                mat_stack.push(temp_mat);
                push_flag = false;
            }
            else
            {
                mat_stack.pop();
                mat_stack.push(temp_mat);
            }
        }

        else if(line == "push")
        {
            //cout<<"push command"<<endl;
            push_flag = true;
        }

        else if(line == "pop")
        {
            //cout<<"pop command"<<endl;
            mat_stack.pop();
        }

        else if(line == "end")
        {
            //cout<<"end command"<<endl;
            break;
        }

        else if(line == "")
        {
            continue;
        }
        else
        {
            for(int i=0;i<line.size();i++)
            {
                printf("p: %c %d\n",line[i],line[i]);
            }
            cout<<"INVALID COMMAND FOUND!!!"<<endl;
        }
    }

    fclose(stdin);
    //stage 2

    ofstream s1("stage1.txt");

    for(int i=0;i<triangles.size();i++)
    {
        s1 << triangles[i].p1.x<<" "<<triangles[i].p1.y<<" "<<triangles[i].p1.z<<endl;
        s1 << triangles[i].p2.x<<" "<<triangles[i].p2.y<<" "<<triangles[i].p2.z<<endl;
        s1 << triangles[i].p3.x<<" "<<triangles[i].p3.y<<" "<<triangles[i].p3.z<<endl;
        s1<<endl;
    }

    s1.close();

    //cout<<"stage 2: "<<endl;
    Point4D l,u,r,l_norm,r_norm,u_norm;
    Point4D minus_eye(-eye.x,-eye.y,-eye.z,1);

    l = look.resultant(minus_eye);
    l_norm = l.normalize();

    r = l.cross(up);
    r_norm = r.normalize();

    u = r_norm.cross(l_norm);

    vector<double> eye_array;
    for(int i=0;i<mat_size;i++)
    {
        for(int j=0;j<mat_size;j++)
        {
            if(i==j) eye_array.push_back(1);
            else if(j==mat_size-1)
            {
                if(i==0) eye_array.push_back(-eye.x);
                else if(i==1) eye_array.push_back(-eye.y);
                else if(i==2) eye_array.push_back(-eye.z);
            }
            else eye_array.push_back(0);
        }
    }
    Matrix eye_mat(eye_array);
    //eye_mat.print();

    vector<double> rot_array;
    for(int i=0;i<mat_size;i++)
    {
        Point4D temp_w(0,0,0,0);
        Point4D minus_l(-l_norm.x,-l_norm.y,-l_norm.z,1);
        for(int j=0;j<mat_size;j++)
        {
            Point4D temp;
            if(i==0) temp = r_norm;
            else if(i==1) temp = u;
            else if(i==2) temp = minus_l;
            else if(i==3) temp = temp_w;

            if(j==0) rot_array.push_back(temp.x);
            else if(j==1) rot_array.push_back(temp.y);
            else if(j==2) rot_array.push_back(temp.z);
            else if(j==3) rot_array.push_back(1-temp.w);
        }
    }
    Matrix cammera(rot_array);
    //cammera.print();

    Matrix view_matrix = cammera.multiply(eye_mat);
    //view_matrix.print();


    for(int i=0;i<triangles.size();i++)
    {
        Point4D p1,p2,p3;
        p1 = view_matrix.multiply_point(triangles[i].p1);
        p2 = view_matrix.multiply_point(triangles[i].p2);
        p3 = view_matrix.multiply_point(triangles[i].p3);

        Triangle new_triangle(p1,p2,p3);
        triangles[i] = new_triangle;
        //new_triangle.print();
    }

    ofstream s2("stage2.txt");

    for(int i=0;i<triangles.size();i++)
    {
        s2 << triangles[i].p1.x<<" "<<triangles[i].p1.y<<" "<<triangles[i].p1.z<<endl;
        s2 << triangles[i].p2.x<<" "<<triangles[i].p2.y<<" "<<triangles[i].p2.z<<endl;
        s2 << triangles[i].p3.x<<" "<<triangles[i].p3.y<<" "<<triangles[i].p3.z<<endl;
        s2<<endl;
    }

    s2.close();


    //stage 3
    //cout<<"stage 3: "<<endl;

    double fovx,t_val,r_val,fovy,ar,near,far;
    fovy = perspective.x;
    ar = perspective.y;
    near = perspective.z;
    far = perspective.w;

    fovx = fovy * ar;
    t_val = near * tan((pi*(fovy*1.0/2.0))/180.0);
    r_val = near * tan((pi*(fovx*1.0/2.0))/180.0);

    vector<double> proj_arr;
    for(int i=0;i<16;i++) proj_arr.push_back(0);
    proj_arr[0] = (near*1.0)/(r_val*1.0);
    proj_arr[5] = (near*1.0)/(t_val*1.0);
    proj_arr[10] = -(far+near*1.0)/(far-near*1.0);
    proj_arr[11] = -(2*far*near*1.0)/(far-near*1.0);
    proj_arr[14] = -1;

    Matrix proj_mat(proj_arr);
    //proj_mat.print();

    for(int i=0;i<triangles.size();i++)
    {
        Point4D p1,p2,p3;
        p1 = proj_mat.multiply_point(triangles[i].p1);
        p2 = proj_mat.multiply_point(triangles[i].p2);
        p3 = proj_mat.multiply_point(triangles[i].p3);

        Triangle new_triangle(p1,p2,p3);
        triangles[i] = new_triangle;
        //new_triangle.print();
    }

    ofstream s3("stage3.txt");

    for(int i=0;i<triangles.size();i++)
    {
        s3 << triangles[i].p1.x<<" "<<triangles[i].p1.y<<" "<<triangles[i].p1.z<<endl;
        s3 << triangles[i].p2.x<<" "<<triangles[i].p2.y<<" "<<triangles[i].p2.z<<endl;
        s3 << triangles[i].p3.x<<" "<<triangles[i].p3.y<<" "<<triangles[i].p3.z<<endl;
        s3<<endl;
    }

    s3.close();

    //stage 4
    string inp;
    ifstream fromFile("./test/1/config.txt");

    getline(fromFile,inp);
    stringstream ss1(inp);
    ss1 >> screen_width;
    ss1 >> screen_height;

    getline(fromFile,inp);
    stringstream ss2(inp);
    ss2 >> left_lim_x;

    getline(fromFile,inp);
    stringstream ss3(inp);
    ss3 >> bottom_lim_y;

    getline(fromFile,inp);
    stringstream ss4(inp);
    ss4 >> front_z;
    ss4 >> rear_z;

    /*cout<<screen_width<<endl;
    cout<<screen_height<<endl;

    cout<<left_lim_x<<endl;
    cout<<bottom_lim_y<<endl;

    cout<<front_z<<endl;
    cout<<rear_z<<endl;*/

    right_lim_x = -left_lim_x;
    top_lim_y = -bottom_lim_y;

    dx = ((right_lim_x-left_lim_x)*1.0)/(screen_width*1.0);
    dy = ((top_lim_y-bottom_lim_y)*1.0)/(screen_height*1.0);

    top_y = top_lim_y - (dy/2.0);
    bottom_y = bottom_lim_y + (dy/2.0);

    left_x = left_lim_x + (dx/2.0);
    right_x = right_lim_x - (dx/2.0);

    Color black(255,255,255);
    vector<vector<double>> z_buffer(screen_width,vector<double>(screen_height,rear_z));
    vector<vector<Color>> frame_buffer(screen_width,vector<Color>(screen_height,black));

    bitmap_image image(screen_width,screen_height);
    for(int i=0;i<triangles.size();i++)
    //for(int i=0;i<1;i++) //for debug only
    {
        //cout<<"=================New Triangle "<<i+1<<"=========="<<endl;
        Point4D point1,point2,point3;
        point1 = triangles[i].p1;
        point2 = triangles[i].p2;
        point3 = triangles[i].p3;

        double top_point = max(max(point1.y,point2.y),point3.y);
        double bottom_point = min(min(point1.y,point2.y),point3.y);
        double left_point = min(min(point1.x,point2.x),point3.x);
        double right_point = max(max(point1.x,point2.x),point3.x);

        double top_scanline = min(top_y,top_point);
        double bottom_scanline = max(bottom_y,bottom_point);

        int y_start = round(((top_y-top_scanline)*1.0)/dy);
        int y_end = round(((top_y-bottom_scanline)*1.0)/dy);  //eta thik ache topy-bottom_scanline hobe
        //cout<<y_start<< "  "<<y_end<<endl;

        for(int j=y_start;j<y_end;j++)
        {
            double scanline_y = top_y - (j*dy);
            pair<double,double> intersect12,intersect13,intersect23;
            double intx,intz;
            double left_scanline,right_scanline;
            int x_start,x_end;
            bool parallel_found = false;
            bool point_found = false;

            if(parallel_found)
            {
                continue;
            }
            else if(point1.isParallel(point2))
            {
                //do later
                if(abs(scanline_y-point1.y) < (dy/2))
                {
                    parallel_found = true;
//                    cout<<j<<endl;
//                    cout<<"here"<<endl;
                    intersect13.first = point1.x;
                    intersect23.first = point2.x;
                    intersect13.second = point1.z;
                    intersect23.second = point2.z;

                    intersect12.first = screen_width*4;
                }
                else
                {
                    intersect12.first= screen_width*2;
                }
            }
            else if(abs(point1.y-scanline_y)<(dy/2))
            {
                //cout<<"p: "<<j<<endl;
                intersect12.first = point1.x;
                intersect12.second = point1.z;
            }

            else
            {
//                cout<<"d_y"<<point1.y-point2.y<<endl;
                intx = (((scanline_y-point1.y)*(point1.x-point2.x)*1.0)/(point1.y-point2.y)) + point1.x;
                intz = (((scanline_y-point1.y)*(point1.z-point2.z)*1.0)/(point1.y-point2.y)) + point1.z;
                intersect12.first=intx;
                intersect12.second=intz;
            }

            if(parallel_found)
            {
                continue;
            }
            else if(point1.isParallel(point3))
            {
                if(abs(scanline_y-point1.y) < (dy/2))
                {
                    parallel_found = true;
//                    cout<<j<<endl;
//                    cout<<"here2"<<endl;
                    intersect12.first = point1.x;
                    intersect23.first = point3.x;
                    intersect12.second = point1.z;
                    intersect23.second = point3.z;

                    intersect13.first = screen_width*4;
                }
                else
                {
                    intersect13.first= screen_width*2;
                }
            }
            else if(abs(point3.y-scanline_y)<(dy/2))
            {
                //cout<<"p: "<<j<<endl;
                intersect13.first = point3.x;
                intersect13.second = point3.z;
            }
            else
            {
                //cout<<"d_y"<<point1.y-point3.y<<endl;
                intx = (((scanline_y-point1.y)*(point1.x-point3.x)*1.0)/(point1.y-point3.y)) + point1.x;
                intz = (((scanline_y-point1.y)*(point1.z-point3.z)*1.0)/(point1.y-point3.y)) + point1.z;
                intersect13.first=intx;
                intersect13.second=intz;
            }

            if(parallel_found)
            {
                continue;
            }
            else if(point2.isParallel(point3))
            {
                if(abs(scanline_y-point2.y) < (dy/2))
                {
                    parallel_found = true;
//                    cout<<j<<endl;
//                    cout<<"here3"<<endl;
                    intersect12.first = point2.x;
                    intersect13.first = point3.x;
                    intersect12.second = point2.z;
                    intersect13.second = point3.z;

                    intersect23.first = screen_width*4;
                }
                else
                {
                    intersect23.first= screen_width*2;
                }
            }
            else if(abs(point2.y-scanline_y)<(dy/2))
            {
                //cout<<"p: "<<j<<endl;
                intersect23.first = point2.x;
                intersect23.second = point2.z;
            }
            else
            {
  //              cout<<"d_y"<<point2.y-point3.y<<endl;
                intx = (((scanline_y-point2.y)*(point2.x-point3.x)*1.0)/(point2.y-point3.y)) + point2.x;
                intz = (((scanline_y-point2.y)*(point2.z-point3.z)*1.0)/(point2.y-point3.y)) + point2.z;
                intersect23.first=intx;
                intersect23.second=intz;
            }
            parallel_found = false;

            vector<pair<double,double>>intersects;
            if(intersect12.first<=right_point && intersect12.first>=left_point)
            {
                intersects.push_back(intersect12);
            }
            if(intersect13.first<=right_point && intersect13.first>=left_point)
            {
                intersects.push_back(intersect13);
            }
            if(intersect23.first<=right_point && intersect23.first>=left_point)
            {
                intersects.push_back(intersect23);
            }

            if(intersects.size() == 3)
            {
                if(abs(intersect12.first-intersect13.first)<(dx/2))
                {
                    intersects.clear();
                    intersects.push_back(intersect12);
                    intersects.push_back(intersect23);
                }
                else if(abs(intersect12.first-intersect23.first)<(dx/2))
                {
                    intersects.clear();
                    intersects.push_back(intersect23);
                    intersects.push_back(intersect13);
                }
                else if(abs(intersect13.first-intersect23.first)<(dx/2))
                {
                    intersects.clear();
                    intersects.push_back(intersect13);
                    intersects.push_back(intersect12);
                }
            }

            /*if(intersects.size() > 2)
            {
                cout<<"error"<<endl;
            }*/

            if(intersects.size() > 1)
            {
                left_scanline = max(min(intersects[0].first,intersects[1].first),left_lim_x);
                right_scanline = min(max(intersects[0].first,intersects[1].first),right_lim_x);
            }
            else if(intersects.size() == 1)
            {
                left_scanline = max(intersects[0].first,left_lim_x);
                right_scanline = min(intersects[0].first,right_lim_x);
            }
            else
            {
                continue;
            }
//            cout<<j<<" here"<<endl;


            x_start = round(((left_scanline-left_x)*1.0)/dx);
            x_end = round(((right_scanline-left_x)*1.0)/dx);

            double zp;
            if(intersects[0].first < intersects[1].first) zp = intersects[0].second;
            else zp = intersects[1].second;

            double diff_x = intersects[0].first - intersects[1].first;
            double diff_z = intersects[0].second - intersects[1].second;

            for(int k=x_start;k<=x_end;k++)
            {
                zp = zp + ((diff_z*1.0)/(diff_x*1.0)) * dx;
                if(zp < z_buffer[j][k] && zp >=  front_z)
                {
                    z_buffer[j][k] = zp;
                    Color cr = triangles[i].c;
                    image.set_pixel(k,j,cr.r,cr.g,cr.b);
                }
            }
        }
    }
    ofstream s4("z_buffer.txt");

    for(int i=0;i<screen_width;i++)
    {
        for(int j=0;j<screen_height;j++)
        {
            if(z_buffer[i][j] < rear_z)
            {
                s4<<z_buffer[i][j]<<"\t";
            }
        }
        s4<<endl;
    }
    s4.close();
    image.save_image("out.bmp");

    return 0;
}

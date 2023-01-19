#include<iostream>
#include<fstream>
#include<math.h>
#include<vector>
#include<omp.h>
#include<array>
#include<algorithm>

class domain{
        
    public:
        double D,D_dowel,LX,LY,LZ;
        int npar, layer;
        std::vector<std::array<double,3>> pos;
    
    public:
        domain(double D, double D_dowel, double LX, double LY, double LZ, int npar, int layer){
            this->D = D;
            this->D_dowel = D_dowel;
            this->LX = LX;
            this->LY = LY;
            this->LZ = LZ;
            this->npar = npar;
            this->layer = layer;
        }

        double n_force(const double x, const double rlimit){
            double cof = 1000/0.05;
            if(x<rlimit-0.05){
                return 1000;
            }
            else if(x<rlimit){
                return cof*(rlimit-x);
            }
            else{
                return 0;
            }
        }

        double minimal_dist(const std::array<double,3> &p1, const std::array<double,3> &p2){
            double dx = abs(p1[0]-p2[0]);
            double dy = abs(p1[1]-p2[1]);
            double dz = abs(p1[2]-p2[2]);
            if(dx>LX/2){
                dx = LX-dx;
            }
            if(dy>LY/2){
                dy = LY-dy;
            }
            if(dz>LZ/2){
                dz = LZ-dz;
            }
            return sqrt(dx*dx+dy*dy+dz*dz);
        }

        std::array<double,3> normal_direction(const std::array<double,3> &p1, const std::array<double,3> &p2){
            double dx = p1[0]-p2[0];
            double dy = p1[1]-p2[1];
            double dz = p1[2]-p2[2];
            if(dx>LX/2){
                dx = dx-LX;
            }
            else if(dx<-LX/2){
                dx = dx+LX;
            }
            if(dy>LY/2){
                dy = dy-LY;
            }
            else if(dy<-LY/2){
                dy = dy+LY;
            }
            if(dz>LZ/2){
                dz = dz-LZ;
            }
            else if(dz<-LZ/2){
                dz = dz+LZ;
            }
            double norm = sqrt(dx*dx+dy*dy+dz*dz);
            return {dx/norm, dy/norm, dz/norm};
        }

        bool check_overlap(const std::vector<std::array<double,3>> &pars, const std::array<double,3> &p, const double rlimit, const int count){
            if (count==0){
                return false;
            }
            for(auto &p1:pars){
                double dist = minimal_dist(p1,p);
                if(dist<rlimit){
                    return true;
                }
            }
            return false;
        }

        std::vector<std::array<double,3>> create_layer_fixz(const int &n, const double &zbase){
            std::vector<std::array<double,3>> pos(n,std::array<double,3>{0,0,0});
            int count = 0;
            while(count<n){
                bool to_generate = true;
                std::array<double,3> p = {LX*rand()/RAND_MAX, LY*rand()/RAND_MAX, zbase};
                while(to_generate){
                    p = {LX*rand()/RAND_MAX, LY*rand()/RAND_MAX, zbase};
                    to_generate = check_overlap(pos,p,0.8,count);
                }
                pos[count] = p;
                count++;
            }

            std::cout<<"layer "<<layer<<" has "<<count<<" particles and has been generated"<<std::endl;
            return pos;
        }

        inline std::array<double,3> scalar_multiply(const std::array<double,3> &p, const double s){
            return {p[0]*s, p[1]*s, p[2]*s};
        }

        inline std::array<double,3> array_minus(const std::array<double,3> &p1, const std::array<double,3> &p2){
            return {p1[0]-p2[0], p1[1]-p2[1], p1[2]-p2[2]};
        }

        inline std::array<double,3> array_plus(const std::array<double,3> &p1, const std::array<double,3> &p2){
            return {p1[0]+p2[0], p1[1]+p2[1], p1[2]+p2[2]};
        }

        inline double array_dot(const std::array<double,3> &p1, const std::array<double,3> &p2){
            return p1[0]*p2[0]+p1[1]*p2[1]+p1[2]*p2[2];
        }

        void move_par_2D(std::vector<std::array<double,3>> &pos, int fixed, double dt, double rlimit){
            int n = pos.size();
            int loop_indice = 0;
            double min_dist=0.0;
            // std::vector<std::array<double,3>> vel(n,std::array<double,3>{0,0,0});
            // std::vector<std::array<double,3>> omega(n,std::array<double,3>{0,0,0});
            // std::vector<std::array<double,3>> force(n,std::array<double,3>{0,0,0});
            // std::vector<std::array<double,3>> torque(n,std::array<double,3>{0,0,0});
            // double min_dist_single[n]={1.0};

            while(min_dist<rlimit){
            //use openmp to parallelize the loop
            std::vector<std::array<double,3>> vel(n,std::array<double,3>{0,0,0});
            std::vector<std::array<double,3>> omega(n,std::array<double,3>{0,0,0});
            std::vector<std::array<double,3>> force(n,std::array<double,3>{0,0,0});
            std::vector<std::array<double,3>> torque(n,std::array<double,3>{0,0,0});
            std::vector<double> min_dist_single(n,1.0);

            #pragma omp parallel for schedule(dynamic,20) 
            for(int i=0;i<n;++i){
                if(i<fixed) continue;

                for (int j=0;j<n;++j){
                    if(i==j) continue;
                    auto nd = normal_direction(pos[j],pos[i]);
                    double dist = minimal_dist(pos[i],pos[j]);
                    double gapp = rlimit-dist;

                    if(dist<min_dist_single[i]){
                        min_dist_single[i] = dist;
                        //std::cout<<i<<" "<<j<<" "<<dist<<pos[i][0]<<" "<<pos[i][1]<<" "<<pos[i][2]<<" "<<pos[j][0]<<" "<<pos[j][1]<<" "<<pos[j][2]<<std::endl;
                    }
                    
                    //std::cout<<i<<" "<<j<<" "<<dist<<" "<<gapp<<std::endl;
                    if (gapp>0)
                    {
                        double rcon1 = rlimit/2 - gapp/2;

                        std::array<double,3> vdiff = {0.0,0.0,0.0};

                        vdiff[0] = vel[i][0]-vel[j][0] + rcon1*(omega[i][1]*nd[2]-omega[i][2]*nd[1])+rcon1*(omega[j][1]*nd[2]-omega[j][2]*nd[1]);
                        vdiff[1] = vel[i][1]-vel[j][1] + rcon1*(omega[i][2]*nd[0]-omega[i][0]*nd[2])+rcon1*(omega[j][2]*nd[0]-omega[j][0]*nd[2]);
                        vdiff[2] = vel[i][2]-vel[j][2] + rcon1*(omega[i][0]*nd[1]-omega[i][1]*nd[0])+rcon1*(omega[j][0]*nd[1]-omega[j][1]*nd[0]);

                        auto vdiff_norm = scalar_multiply(nd,nd[0]*vdiff[0]+nd[1]*vdiff[1]+nd[2]*vdiff[2]);
                        auto vdiff_tan  = array_minus(vdiff,vdiff_norm);
                        auto vtt =std::sqrt(vdiff_tan[0]*vdiff_tan[0]+vdiff_tan[1]*vdiff_tan[1]+vdiff_tan[2]*vdiff_tan[2]);

                        std::array<double,3> cospt = {0.0,0.0,0.0};
                        if(vtt>0) cospt = scalar_multiply(vdiff_tan,1.0/vtt);
 
                        auto nf = n_force(dist,rlimit);
                        force[i][0] += -1.0*nd[0]*nf*0.004 - vdiff_norm[0]*0.01;
                        force[i][1] += -1.0*nd[1]*nf*0.004 - vdiff_norm[1]*0.01;
                        force[i][0] += -vdiff_tan[0]*0.01*cospt[0];
                        force[i][1] += -vdiff_tan[1]*0.01*cospt[1];

                        force[i][2] = 0.0;

                        torque[i][0] += nd[2]*vdiff_tan[1] - nd[1]*vdiff_tan[2];
                        torque[i][1] += nd[0]*vdiff_tan[2] - nd[2]*vdiff_tan[0];
                        torque[i][2] += nd[1]*vdiff_tan[0] - nd[0]*vdiff_tan[1];

                    }
                }
            }

            min_dist = *std::min_element(min_dist_single.begin(),min_dist_single.end());

            if(min_dist>rlimit) {
                std::cout<<"loop:"<<loop_indice<<"with min_dist:"<<min_dist<<std::endl;
                break;
            }

            if(loop_indice>1000){
                std::cout<<"too many loops "<<std::endl;
                break;
            }
            loop_indice++;

            #pragma omp parallel for schedule(dynamic,20)
            for(int i=0;i<n;++i){
                if(i<fixed) continue;
                vel[i] = array_plus(vel[i],scalar_multiply(force[i],dt));
                pos[i] = array_plus(pos[i],scalar_multiply(vel[i],dt));
                omega[i] = array_plus(omega[i],scalar_multiply(torque[i],0.001*dt));
                //if(i==100&&loop_indice%100==0) std::cout<<i<<" "<<pos[i][0]<<" "<<pos[i][1]<<" "<<vel[i][0]<<" "<<vel[i][1]<<" "<<scalar_multiply(force[i],dt)[0]<<" "<<force[i][0]<<std::endl;
            }

            //for periodic boundary condition
            #pragma omp parallel for schedule(dynamic,20)
            for(int i=0;i<n;++i){
                if(i<fixed) continue;
                if(pos[i][0]>LX) pos[i][0] -= LX;
                if(pos[i][0]<0) pos[i][0] += LX;
                if(pos[i][1]>LY) pos[i][1] -= LY;
                if(pos[i][1]<0) pos[i][1] += LY;
            }

        }
        std::cout<<"target:"<<rlimit<<"loop"<<loop_indice<<"with min_dist:"<<min_dist<<std::endl;
        }

        void write_to_file(std::string filename){
            std::ofstream file;
            file.open(filename);
            for(int i=0;i<npar;++i){
                file<<pos[i][0]<<" "<<pos[i][1]<<" "<<pos[i][2]<<std::endl;
            }
            file.close();
        }

        void add_to_pos(std::vector<std::array<double,3>> add){
            for(auto &a:add){
                pos.push_back(a);
            }
        }
};

int main(){

    domain dom(1,1.15,60,30,10,191*90,10);

    //auto pos = dom.create_layer_fixz(191*9,0.5);

    for(int l=0;l<dom.layer;++l){
        auto pos = dom.create_layer_fixz(dom.npar/dom.layer,0.5+l*1.0);
        
        for(int i=0;i<21;++i){
            dom.move_par_2D(pos,0,0.02,0.81+i*0.01);
        }

        dom.add_to_pos(pos);
    }

    dom.write_to_file("pos.txt");

    return 0;

    
}

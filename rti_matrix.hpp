// author: Jungwook Shin @ mgh (jshin13@mgh.harvard.edu)

#ifndef RTI_MATRIX_H
#define RTI_MATRIX_H

#include <iostream>
#include <cmath>
#include <array>

#include "rti_vec.hpp"

namespace rti{

//interpolation function::ICRU/NIST:range energy
//this is a copy of HepRotation.
//public member has no '_' at the end of its name

template<typename T>
class mat3x3 {
public:
    //angle
    T x;
    T y;
    T z;
    //matrix element
    T xx;
    T xy;
    T xz;
    T yx;
    T yy;
    T yz;
    T zx;
    T zy;
    T zz;
    
    CUDA_HOST_DEVICE
    mat3x3(): 
    x(0),y(0),z(0),
    xx(1.0),xy(0),xz(0),
    yx(0),yy(1.0),yz(0),
    zx(0),zy(0),zz(1.0)
    {;}

    CUDA_HOST_DEVICE
    mat3x3(T xx, T xy, T xz,T yx, T yy, T yz,T zx, T zy, T zz ): 
    x(0),y(0),z(0),
    xx(xx),xy(xy),xz(xz),
    yx(yx),yy(yy),yz(yz),
    zx(zx),zy(zy),zz(zz)
    {;}

    CUDA_HOST_DEVICE
    mat3x3(const mat3x3& ref){
        x  = ref.x ;
        y  = ref.y;
        z  = ref.z;
        xx = ref.xx;
        xy = ref.xy;
        xz = ref.xz;
        yx = ref.yx;
        yy = ref.yy;
        yz = ref.yz;
        zx = ref.zx;
        zy = ref.zy;
        zz = ref.zz;
    }

    CUDA_HOST_DEVICE
    mat3x3(T a, T b, T c):
        x(a),y(b),z(c),
        xx(1.0),xy(0),xz(0),
        yx(0),yy(1.0),yz(0),
        zx(0),zy(0),zz(1.0){
        if (x != 0) this->rotate_x(x);
        if (y != 0) this->rotate_y(y);
        if (z != 0) this->rotate_z(z);
    }

    CUDA_HOST_DEVICE
    mat3x3(std::array<T,3>& abc): 
        x(abc[0]),y(abc[1]),z(abc[2]),
        xx(1.0),xy(0),xz(0),
        yx(0),yy(1.0),yz(0),
        zx(0),zy(0),zz(1.0)
    {
        if (x != 0) this->rotate_x(x);
        if (y != 0) this->rotate_y(y);
        if (z != 0) this->rotate_z(z);
    }

    CUDA_HOST_DEVICE
    ~mat3x3(){;}

    CUDA_HOST_DEVICE
    mat3x3& 
    rotate_x(T a){
        x = a;
        #if defined(__CUDACC__)
        T c1 = cosf(x);
        T s1 = sinf(x);
        #else
        T c1 = std::cos(x);
        T s1 = std::sin(x);
        #endif
        T x1 = yx, y1 = yy, z1 = yz; 
        yx = c1*x1 - s1*zx;
        yy = c1*y1 - s1*zy;
        yz = c1*z1 - s1*zz;
        zx = s1*x1 + c1*zx;
        zy = s1*y1 + c1*zy;
        zz = s1*z1 + c1*zz;
        return *this;
    }

    CUDA_HOST_DEVICE
    mat3x3&
    rotate_y(T a){
        y = a;
        #if defined(__CUDACC__)
        T c1 = cosf(y);
        T s1 = sinf(y);
        #else
        T c1 = std::cos(y);
        T s1 = std::sin(y);
        #endif
        T x1 = zx, y1 = zy, z1 = zz; 
        zx = c1*x1 - s1*xx;
        zy = c1*y1 - s1*xy;
        zz = c1*z1 - s1*xz;
        xx = s1*x1 + c1*xx;
        xy = s1*y1 + c1*xy;
        xz = s1*z1 + c1*xz;
        return *this;
    }
    
    CUDA_HOST_DEVICE
    mat3x3&
    rotate_z(T a){
        z = a;
        #if defined(__CUDACC__)
        T c1 = cosf(z);
        T s1 = sinf(z);
        #else
        T c1 = std::cos(z);
        T s1 = std::sin(z);
        #endif
        T x1 = xx, y1 = xy, z1 = xz; 
        xx = c1*x1 - s1*yx;
        xy = c1*y1 - s1*yy;
        xz = c1*z1 - s1*yz;
        yx = s1*x1 + c1*yx;
        yy = s1*y1 + c1*yy;
        yz = s1*z1 + c1*yz;
        return *this;
    }

    CUDA_HOST_DEVICE
    std::array<T,3>
    operator * (const std::array<T,3>& r) const {
        return std::array<T,3>({ 
            xx * r[0] + xy * r[1] + xz * r[2],
            yx * r[0] + yy * r[1] + yz * r[2],
            zx * r[0] + zy * r[1] + zz * r[2]});
    }

    CUDA_HOST_DEVICE
    rti::vec3<T>
    operator * (const rti::vec3<T>& r) const {
        return rti::vec3<T>(
            xx * r.x + xy * r.y + xz * r.z,
            yx * r.x + yy * r.y + yz * r.z,
            zx * r.x + zy * r.y + zz * r.z);
    }

    CUDA_HOST_DEVICE
    mat3x3 inverse() const {        
        return rti::mat3x3<T>(xx, yx, zx, xy, yy, zy, xz, yz, zz );
    }

    /*
    mat3x3<T> 
    operator* (const mat3x3<T>& r){
        return mat3x3<T>(
            xx*r.xx + xy*r.yx + xz*r.zx,
            xx*r.xy + xy*r.yy + xz*r.zy,
            xx*r.xz + xy*r.yz + xz*r.zz,
            yx*r.xx + yy*r.yx + yz*r.zx,
            yx*r.xy + yy*r.yy + yz*r.zy,
            yx*r.xz + yy*r.yz + yz*r.zz,
            zx*r.xx + zy*r.yx + zz*r.zx,
            zx*r.xy + zy*r.yy + zz*r.zy,
            zx*r.xz + zy*r.yz + zz*r.zz );
    }*/

    CUDA_HOST_DEVICE
    mat3x3<T>
    operator* (const mat3x3<T>& r) const {
        return mat3x3<T>(
            xx*r.xx + xy*r.yx + xz*r.zx,
            xx*r.xy + xy*r.yy + xz*r.zy,
            xx*r.xz + xy*r.yz + xz*r.zz,
            yx*r.xx + yy*r.yx + yz*r.zx,
            yx*r.xy + yy*r.yy + yz*r.zy,
            yx*r.xz + yy*r.yz + yz*r.zz,
            zx*r.xx + zy*r.yx + zz*r.zx,
            zx*r.xy + zy*r.yy + zz*r.zy,
            zx*r.xz + zy*r.yz + zz*r.zz );
    }

    CUDA_HOST_DEVICE
    void dump() const {
    #if defined(__CUDACC__)
        printf("xx,xy,xz %f, %f, %f\n", xx,xy, xz);
        printf("yx,yy,yz %f, %f, %f\n", yx,yy, yz);
        printf("zx,zy,zz %f, %f, %f\n", zx,zy, zz);
    #else
        std::cout<<"xx,xy,xz "<< xx <<" " << xy <<" " << xz << std::endl;
        std::cout<<"yx,yy,yz "<< yx <<" " << yy <<" " << yz << std::endl;
        std::cout<<"zx,zy,zz "<< zx <<" " << zy <<" " << zz << std::endl;
    #endif
    }
};


template<typename T>
class mat4x4 {
public:
    //matrix element
    T xx;
    T xy;
    T xz;
    T xs;
    T yx;
    T yy;
    T yz;
    T ys;
    T zx;
    T zy;
    T zz;
    T zs;
    T sx;
    T sy;
    T sz;
    T ss;
    
    CUDA_HOST_DEVICE
    mat4x4(): 
    xx(1.0),xy(0),xz(0),xs(0),
    yx(0),yy(1.0),yz(0),ys(0),
    zx(0),zy(0),zz(1.0),zs(0),
    sx(0),sy(0),sz(0.0),ss(1.0)
    {;}

    CUDA_HOST_DEVICE
    mat4x4(
        T xx, T xy, T xz, T xs,
        T yx, T yy, T yz, T ys,
        T zx, T zy, T zz, T zs,
        T sx, T sy, T sz, T ss
    ): 
    xx(xx),xy(xy),xz(xz),xs(xs),
    yx(yx),yy(yy),yz(yz),ys(ys),
    zx(zx),zy(zy),zz(zz),zs(zs),
    sx(sx),sy(sy),sz(sz),ss(ss)
    {;}

    CUDA_HOST_DEVICE
    mat4x4(const mat4x4& ref){
        xx = ref.xx;
        xy = ref.xy;
        xz = ref.xz;
        xs = ref.xs;

        yx = ref.yx;
        yy = ref.yy;
        yz = ref.yz;
        ys = ref.ys;
        
        zx = ref.zx;
        zy = ref.zy;
        zz = ref.zz;
        zs = ref.zs;

        sx = ref.sx;
        sy = ref.sy;
        sz = ref.sz;
        ss = ref.ss;
    }

    CUDA_HOST_DEVICE
    mat4x4(const T* a)
    {
        xx = a[0]; xy = a[1]; xz = a[2]; xs = a[3];
        yx = a[4]; yy = a[5]; yz = a[6]; ys = a[7];
        zx = a[8]; zy = a[9]; zz = a[10]; zs = a[11];
        sx = a[12]; sy = a[13]; sz = a[14];ss = a[15];
    }

    CUDA_HOST_DEVICE
    ~mat4x4(){;}

    CUDA_HOST_DEVICE
    std::array<T,4>
    operator * (const std::array<T,4>& r) const {
        return std::array<T,4>({ 
            xx * r[0] + xy * r[1] + xz * r[2] + xs*r[3],
            yx * r[0] + yy * r[1] + yz * r[2] + ys*r[3],
            zx * r[0] + zy * r[1] + zz * r[2] + zs*r[3],
            sx * r[0] + sy * r[1] + sz * r[2] + ss*r[3]});
    }

    CUDA_HOST_DEVICE
    rti::vec4<T>
    operator * (const rti::vec4<T>& r) const {
        return rti::vec4<T>(
            xx * r.x + xy * r.y + xz * r.z + xs * r.s,
            yx * r.x + yy * r.y + yz * r.z + ys * r.s,
            zx * r.x + zy * r.y + zz * r.z + zs * r.s,
            sx * r.x + sy * r.y + sz * r.z + ss * r.s);
    }

    /*
    mat4x4<T> 
    operator* (const mat3x3<T>& r){
        return mat3x3<T>(
            xx*r.xx + xy*r.yx + xz*r.zx,
            xx*r.xy + xy*r.yy + xz*r.zy,
            xx*r.xz + xy*r.yz + xz*r.zz,
            yx*r.xx + yy*r.yx + yz*r.zx,
            yx*r.xy + yy*r.yy + yz*r.zy,
            yx*r.xz + yy*r.yz + yz*r.zz,
            zx*r.xx + zy*r.yx + zz*r.zx,
            zx*r.xy + zy*r.yy + zz*r.zy,
            zx*r.xz + zy*r.yz + zz*r.zz );
    }*/

    CUDA_HOST_DEVICE
    mat4x4<T>
    operator* (const mat4x4<T>& r) const {
        return mat4x4<T>(
            xx*r.xx + xy*r.yx + xz*r.zx + xs*r.sx,
            xx*r.xy + xy*r.yy + xz*r.zy + xs*r.sy,
            xx*r.xz + xy*r.yz + xz*r.zz + xs*r.sz,
            xx*r.xs + xy*r.ys + xz*r.zs + xs*r.ss,
            
            yx*r.xx + yy*r.yx + yz*r.zx + ys*r.sx,
            yx*r.xy + yy*r.yy + yz*r.zy + ys*r.sy,
            yx*r.xz + yy*r.yz + yz*r.zz + ys*r.sz,
            yx*r.xs + yy*r.ys + yz*r.zs + ys*r.ss,
            
            zx*r.xx + zy*r.yx + zz*r.zx + zs*r.sx,
            zx*r.xy + zy*r.yy + zz*r.zy + zs*r.sy,
            zx*r.xz + zy*r.yz + zz*r.zz + zs*r.sz,
            zx*r.xs + zy*r.ys + zz*r.zs + zs*r.ss,
            
            sx*r.xx + sy*r.yx + sz*r.zx + ss*r.sx,
            sx*r.xy + sy*r.yy + sz*r.zy + ss*r.sy,
            sx*r.xz + sy*r.yz + sz*r.zz + ss*r.sz,
            sx*r.xs + sy*r.ys + sz*r.zs + ss*r.ss
            );
    }

    /*
    CUDA_HOST_DEVICE
    mat4x4 inverse() const {        
        T det = xx*(yy*zz*ss + yx*zs*sy + )
        //return rti::mat3x3<T>(xx, yx, zx, xy, yy, zy, xz, yz, zz );
    }
    */

    CUDA_HOST_DEVICE
    void dump() const {
    #if defined(__CUDACC__)
        printf("xx,xy,xz,xs %f, %f, %f, %f\n", xx,xy, xz, xs);
        printf("yx,yy,yz,ys %f, %f, %f, %f\n", yx,yy, yz, ys);
        printf("zx,zy,zz,zs %f, %f, %f, %f\n", zx,zy, zz, zs);
        printf("sx,sy,sz,ss %f, %f, %f, %f\n", sx,sy, sz, ss);
    #else
        std::cout<<"xx,xy,xz,xs "<< xx <<" " << xy <<" " << xz << " " << xs << std::endl;
        std::cout<<"yx,yy,yz,ys "<< yx <<" " << yy <<" " << yz << " " << ys << std::endl;
        std::cout<<"zx,zy,zz,zs "<< zx <<" " << zy <<" " << zz << " " << zs << std::endl;
        std::cout<<"zx,zy,zz,ss "<< sx <<" " << sy <<" " << sz << " " << ss << std::endl;
    #endif
    }
};



}
#endif
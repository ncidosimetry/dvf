// author: Jungwook Shin @ mgh (jshin13@mgh.harvard.edu)

#ifndef RTI_CT_H
#define RTI_CT_H

/*
#include "gdcmDataSet.h"
#include "gdcmDataElement.h"
#include "gdcmVR.h"
#include "gdcmVM.h"
#include "gdcmByteValue.h"
#include "gdcmItem.h"
#include "gdcmGlobal.h"
#include "gdcmDicts.h"
#include "gdcmSequenceOfItems.h"
*/

#include "gdcmScanner.h"
#include "gdcmAttribute.h"
#include "gdcmDirectory.h"
#include "gdcmIPPSorter.h"
#include "gdcmImageReader.h"

#include "rti_rect3d.hpp"
#include "rti_matrix.hpp"

namespace rti{

/**
 * The <code>ct</code> class represents a 3D CT to hold HU per pixel.
 *
 * @tparam R for grid coordinates, float or double.
 */
template<typename R> //T should be HU
class ct : public rect3d<int16_t, R> {
protected:
    std::vector<std::string>   files_            ;  //sorted file in Z acsending order 
    std::map<std::string, std::string> uid2file_ ;  //sop instance UID to id in files_
    char* ct_dir;
    R dx_;
    R dy_;
public: 
    CUDA_HOST
    ct(){;}

    /**
    * Construct a rectlinear grid from array of x/y/z with their size
    * x,y dimension is unchanged but z pixel size may change
    * @param[in] f CT directory.
    */
    CUDA_HOST
    ct(std::string f, bool is_print=false){
        ct_dir = new char[f.length()+1];
        strcpy(ct_dir, f.c_str());

        //http://gdcm.sourceforge.net/html/SortImage_8cxx-example.html#_a5
        gdcm::Directory dir;
        dir.Load(ct_dir, false); //non-recursive
        
        const gdcm::Directory::FilenamesType& all_files = dir.GetFilenames();
    
        gdcm::Scanner scanner ;
        const gdcm::Tag ct_tag = gdcm::Tag(0x08, 0x60);
        scanner.AddTag(ct_tag);
        scanner.Scan(all_files);

        auto ct_files = scanner.GetAllFilenamesFromTagToValue(ct_tag, "CT" );
        
        gdcm::IPPSorter ippsorter;
        ippsorter.SetComputeZSpacing(false);
        ippsorter.Sort(ct_files); //asending along z
        files_ = ippsorter.GetFilenames();

        gdcm::Scanner s;
        s.AddTag( gdcm::Tag(0x0008,0x0018) ); /**< SOP instance UID */
        s.AddTag( gdcm::Tag(0x0020,0x0032) ); /**< Image Position (Patient) */
        s.AddTag( gdcm::Tag(0x0028,0x0010) ); /**< Rows */
        s.AddTag( gdcm::Tag(0x0028,0x0011) ); /**< Columns */
        s.AddTag( gdcm::Tag(0x0028,0x0030) ); /**< Pixel spacing */

        if (!s.Scan( files_ )){
            assert("scan fail.");
        }
    
        size_t nx ;                 /**< columns */
        size_t ny ;                 /**< rows */
        size_t nz = files_.size() ; /**< number of images. */
        rect3d<int16_t, R>::z_ = new R[nz];
        rect3d<int16_t, R>::dim_.z = nz;
        double x0, y0;
        
        for(size_t i = 0 ; i < nz ; ++i){
            gdcm::Scanner::TagToValue const &m0 = s.GetMapping(files_[i].c_str());

            std::string img_position(m0.find(gdcm::Tag(0x0020,0x0032))->second);
            unsigned int deli0 = img_position.find_first_of( '\\' );
            unsigned int deli1 = img_position.find_last_of( '\\' );
            rect3d<int16_t, R>::z_[i] = (R) (std::stod( img_position.substr(deli1+1) ));
    
            /**< We only determine rows, colums, x0, y0, dx, and dy with first image */
            if(i==0){
                x0 = std::stod(img_position.substr(0,deli0));
                y0 = std::stod(img_position.substr(deli0+1, deli1-deli0));

                ny = std::stoi(m0.find(gdcm::Tag(0x0028,0x0010))->second);
                nx = std::stoi(m0.find(gdcm::Tag(0x0028,0x0011))->second); 
                rect3d<int16_t, R>::dim_.x = nx;
                rect3d<int16_t, R>::dim_.y = ny;
                std::string  pixel_spacing(m0.find(gdcm::Tag(0x0028,0x0030))->second) ;
                unsigned int deli = pixel_spacing.find_first_of( '\\' );
                dx_ = std::stod(pixel_spacing.substr(0, deli));
                dy_ = std::stod(pixel_spacing.substr(deli+1));

            }
            
            /**< A map to search file path upon instance UID */
            uid2file_.insert(std::make_pair(m0.find(gdcm::Tag(0x0008,0x0018))->second, files_[i]));
        }

        rect3d<int16_t, R>::x_ = new R[nx];
        for(size_t i = 0 ; i < nx ; ++i){
            rect3d<int16_t, R>::x_[i] = x0 + dx_ * i;
        }
        rect3d<int16_t, R>::y_ = new R[nx];
        for(size_t i = 0 ; i < ny ; ++i){
            rect3d<int16_t, R>::y_[i] = y0 + dy_ * i;
        }
        std::cout<<"CT (nx,ny,nz): (" << nx << ", " << ny <<", " << nz <<")" << std::endl;
        std::cout<<"CT (dx,dy,dz): (" << dx_ << ", " << dy_ <<", " << rect3d<int16_t, R>::z_[1] - rect3d<int16_t, R>::z_[0] <<")" << std::endl;
        std::cout<<"CT (x,y,z): (" << rect3d<int16_t, R>::x_[0] 
                                   << ", " << rect3d<int16_t, R>::y_[0] 
                                   <<", "  << rect3d<int16_t, R>::z_[0] << ")" << std::endl;

    }

    /**
    * Load patient's image to volume
    *
    * @param[in] f CT directory.
    */
    CUDA_HOST
    virtual 
    void
    load_data(){
        size_t nb_voxels_2d = rect3d<int16_t,R>::dim_.x * rect3d<int16_t,R>::dim_.y ;
        size_t nb_voxels_3d = nb_voxels_2d * rect3d<int16_t,R>::dim_.z;

        rect3d<int16_t,R>::data_.resize(nb_voxels_3d);
        float intercept = 0;
        float slope     = 1;

        for(size_t i=0 ; i <  rect3d<int16_t,R>::dim_.z; ++i){
            
            gdcm::ImageReader reader;    
            reader.SetFileName(files_[i].c_str());
            reader.Read();
            const gdcm::Image& img  = reader.GetImage();
            
            //n_x * n_y * bytes = img.GetBufferLength()
            intercept = float(img.GetIntercept());
            slope     = float(img.GetSlope());
            
            gdcm::PixelFormat pixeltype = img.GetPixelFormat();
            
            switch( pixeltype )
            {
            /*
            case gdcm::PixelFormat::INT8:
                {
                std::valarray<int8_t> d(nb_voxels) ;
                img.GetBuffer((char*) &d[0]);
                this->transform_to_float(d, intercept, slope);
                }
                break;
            case gdcm::PixelFormat::UINT8:
                {
                std::valarray<uint8_t> d(nb_voxels) ;
                img.GetBuffer((char*) &d[0]);
                this->transform_to_float(d, intercept, slope);
                }
                break;
            */
            case gdcm::PixelFormat::INT16:
                {
                
                img.GetBuffer((char*) &rect3d<int16_t,R>::data_[i*nb_voxels_2d]);
                
                }
                break;
            /*
            case gdcm::PixelFormat::UINT16:
                {
                std::valarray<uint16_t> d(nb_voxels) ;
                img.GetBuffer((char*) &d[0]);
                this->transform_to_float(d, intercept, slope);
                }
                break;
            case gdcm::PixelFormat::INT32:
                {
                std::valarray<int32_t> d(nb_voxels) ;
                img.GetBuffer((char*) &d[0]);
                this->transform_to_float(d, intercept, slope);
                }
                break;
            case gdcm::PixelFormat::UINT32:
                {
                //Only this Pixelformat was tested.
                std::valarray<int16_t> d(nb_voxels) ;
                img.GetBuffer((char*) &d[0]);
                //img.GetBuffer( (char*) &rect3d<float,R>::data_[0] ); //didn't work
                
                this->transform_to_float(d, intercept, slope);

                }
                break;
            */
            default:
                assert(0);
            }//switch

        }//for

        rect3d<int16_t,R>::data_ = rect3d<int16_t,R>::data_*slope + intercept;
        
    }

    /**
     *  Finding index is simplified cause pixel spacing is same along x
     */
    CUDA_HOST_DEVICE
    inline virtual size_t
    find_c000_x_index(
        const R& x)
    {   
        assert( (x >= rect3d<int16_t,R>::x_[0]) && (x < rect3d<int16_t,R>::x_[rect3d<int16_t,R>::dim_.x-1]) );
        assert( dx_ > 0 );
        return floor(x/dx_);
    }

    /**
     *  Finding index is simplified cause pixel spacing is same along y
     */
    CUDA_HOST_DEVICE
    inline virtual size_t
    find_c000_y_index(
        const R& y)
    {
        assert( (y >= rect3d<int16_t,R>::y_[0]) && (y < rect3d<int16_t,R>::y_[rect3d<int16_t,R>::dim_.y-1]) );
        assert( dy_ > 0 );
        return floor(y/dy_);
    }

    /**
    * A friend function to copy CT's grid information of src to dest
    *
    */
    template<typename R0, typename T1, typename R1>
    friend void clone_ct_structure(ct<R0>& src, rect3d<T1,R1>& dest);
};


}

#endif
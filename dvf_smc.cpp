// author: Jungwook Shin @ mgh (jshin13@mgh.harvard.edu)


#include "../cli_rti.hpp"

#include "rti_rect3d.hpp"
#include "rti_ct.hpp"
#include "rti_dvf.hpp"
#include "rti_rtdose.hpp"

#include <set>

using namespace std;


//command 2.
//
int main(int argc, char** argv){

    rti::cli_dvf cl_opts;
    cl_opts.read(argc, argv);

    /*
    when float as used for pixel_pitch_type, the strange distribution I got
    */
    typedef double pixel_pitch_type;
    //typedef float pixel_pitch_type;

    //ct directory
    std::string dname = cl_opts["--ctdir"][0];
    rti::ct<pixel_pitch_type> myct(dname); 
    myct.load_data();

    //rtdose file
    std::string fname = cl_opts["--dosefile"][0];
    rti::rtdose<pixel_pitch_type> mydose(fname);
    std::cout<<" RTDOSE: -------"<<std::endl;
    mydose.load_data();
    mydose.write_data("output_dose.raw");

    //dvf file
    std::string dvfname = cl_opts["--dvffile"][0];
    rti::dvf<pixel_pitch_type, pixel_pitch_type> mydvf(dvfname);
    std::cout<<" DVF: -------"<<std::endl;
    mydvf.load_data();

    //Clone various grid options for interpolation
    rti::rect3d<float, pixel_pitch_type>  dose0;
    rti::dvf<pixel_pitch_type, pixel_pitch_type>  dvf5to0;
    rti::clone_structure(mydose,   dose0); //
    rti::clone_structure(mydose, dvf5to0); //careful dx_, dy_, dz_ are not copied.

    float zero_dose  = 0.0;
    rti::vec3<pixel_pitch_type> zero_vector(0.0,0.0,0.0);

    //interpolate DVF on fine grid
    rti::interpolate(mydvf, dvf5to0, zero_vector); 
    
    //rti::warp_linear(mydose, dose0, dvf5to0, zero_dose); 
    //this doesn't seem to work...because dx_,dy_,dz_ are not copied..? 
    //dvf5to0.write_data("dvf_copy2.raw");
    //warp with original warping method takes just 2 sec

    rti::warp_linear(mydose, dose0, mydvf, zero_dose);
    dose0.write_data("dose_warp1.raw");
    
    //rti::warp_linear(mydose, dose0, dvf5to0); 
    //interpolated DVF doesn't work for this 
    //dose0.write_data("dose_warp0.raw");
    
    
    return 0;
}
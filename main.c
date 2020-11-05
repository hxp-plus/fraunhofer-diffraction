#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

double slit_width_mm=0.02;            // 单缝宽度单位为mm
double screen_width_mm=0.1;           // 远处接收衍射图像屏幕的宽度，单位为mm
double wavelength_nm=400;             // 光的波长，单位为nm
double y_screen_dm=20;                // 屏幕距离单缝的距离，单位为dm
double integrate_step=0.0000001;      // 积分步长，单位为m
double speed_of_light=300000000.0;    // 光速，单位为m/s

// 计算单缝中一个点x_point，在屏幕上另一点x_screen，时间time时，产生的复振幅
gsl_complex calculateComplexAmplitudeOnePoint(double x_point, double x_screen, double time) {
    double delta_x=fabs(x_point - x_screen);
    double distance=gsl_complex_abs(gsl_complex_rect(delta_x,y_screen_dm*pow(10,-1)));
    // 计算x_point到x_screen有多少波长
    double number_of_wavelength=distance/(wavelength_nm*pow(10,-9));
    // 计算复振幅的角度，需要加上因为时间的流逝，光从单缝射出已经产生的相位
    double amplitude_angle=2*M_PI*((number_of_wavelength)+time*speed_of_light/(wavelength_nm*pow(10,-9)));
    // 返回复振幅，角度为刚刚计算出来的角度，大小为1
    return gsl_complex_polar(1, amplitude_angle);
}

// 计算一个单缝，位置为slit_position[0]到slit_position[1]，到屏幕上x_screen位置，时间为time时，产生的总复振幅
gsl_complex calculateComplexAmplitudeOneSlit(const double slit_position[], double x_screen, double time) {
    gsl_complex total_complex_amplitude=gsl_complex_rect(0,0);
    // 将缝拆分成100个点，然后用上面定义的函数计算复振幅，之后相加，再返回
    for(double i=slit_position[0];i<=slit_position[1];i+=(slit_position[1]-slit_position[0])/1000) {
        total_complex_amplitude =
                gsl_complex_add(calculateComplexAmplitudeOnePoint(i, x_screen, time), total_complex_amplitude);
    }
    return total_complex_amplitude;
}

// 计算双缝光强，两个缝的位置分别是slit_position_a和slit_position_b
gsl_complex calculateComplexAmplitudeTwoSlits(const double slit_position_a[], const double slit_position_b[],
                                              double x_screen, double time) {
    return gsl_complex_add(calculateComplexAmplitudeOneSlit(slit_position_a, x_screen, time),
                           calculateComplexAmplitudeOneSlit(slit_position_b, x_screen, time));
}

int main() {
    clock_t begin = clock();
    printf("程序已启动，开始计时...\n");

    FILE *fp;
    // 将数据保存在/tmp/data.dat这个文件下
    fp=fopen("/tmp/data.dat", "w");
    // 定义单缝的起始位置和终止位置
    double slit_position_a[2] = {-5*slit_width_mm/2.0*pow(10,-3),-3*slit_width_mm/2.0*pow(10,-3)};
    double slit_position_b[2] = {3*slit_width_mm/2.0*pow(10,-3),5*slit_width_mm/2.0*pow(10,-3)};
    // 计算光的一个周期的时间
    double period=wavelength_nm*pow(10,-9)/speed_of_light;
    // 对全屏幕进行积分，写入屏幕上每一点的光强
    for(double y=-screen_width_mm/2.0; y<screen_width_mm/2.0; y+=integrate_step) {
        // 对屏幕上一个确定的点，计算一个周期内不同时间内的光强，再相加，之后计算平均光强，写入数据
        double light_intensity = 0;
        for(double time=0; time<period; time+=integrate_step) {
            light_intensity+=gsl_pow_int(GSL_REAL(
                    calculateComplexAmplitudeTwoSlits(slit_position_a,slit_position_b,y,time)),2);
        }
        fprintf(fp,"%f\t%f\n", y, light_intensity/(period/integrate_step));
    }
    fclose(fp);

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("程序已完成数据处理，总共使用的时间为：%f 秒\n", time_spent);
    printf("正在绘图...\n");
//    system("gnuplot -e \"plot '/tmp/data.dat' with points\" -p");
    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("程序已完成绘图，总共使用的时间为：%f 秒\n", time_spent);
    return 0;
}
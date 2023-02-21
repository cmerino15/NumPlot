#include "NumCpp.hpp"
#include "SupportLib.hpp"
#include "pbPlots.hpp"
#include <cstdlib>
#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
 
using FunctionType = std::function<double(const nc::NdArray<double>&, const nc::NdArray<double>&)>;
void wikipediaExample()
{
    // https://en.wikipedia.org/wiki/Gauss%E2%80%93Newton_algorithm
 
    // In a biology experiment studying the relation between substrate concentration [S] and reaction rate in
    // an enzyme-mediated reaction, the data in the following table were obtained.
    nc::NdArray<double> sMeasured = { 0.038, 0.194, 0.425, 0.626, 1.253, 2.5, 3.74 };
    nc::NdArray<double> rateMeasured = { 0.05, 0.127, 0.094, 0.2122, 0.2729, 0.2665, 0.3317 };
    
    std::vector<double> x = { 0.038, 0.194, 0.425, 0.626, 1.253, 2.5, 3.74 };
    std::vector<double> y =  { 0.05, 0.127, 0.094, 0.2122, 0.2729, 0.2665, 0.3317 };
        
    RGBABitmapImageReference *imageReference = CreateRGBABitmapImageReference();
    double Mx=0, mx =0, my=100, My=0;
        Mx = *max_element(x.begin(), x.end());
        mx = *min_element(x.begin(), x.end());
        // alternatively - careful with choosing my and My
        for( auto i: y){
            My = std::max(My, i);
            my = std::min(my, i);  
        }   
        std::vector<double> xs2 = {mx-1, Mx+1};
        std::vector<double> ys2 = {my-1, My+1};
	ScatterPlotSeries *series = GetDefaultScatterPlotSeriesSettings();
	series->xs = &x;
	series->ys = &y;
	series->linearInterpolation = false;
	series->pointType = toVector(L"dots");
	series->color = CreateRGBColor(.2, .8, .5);

        //calculate linear regression line and replace with min and max values below
        std::vector<double> xplot = {mx, Mx};
        std::vector<double> yplot = {my, My};
        
	ScatterPlotSeries *series2 = GetDefaultScatterPlotSeriesSettings();
	series2->xs = &xplot;
	series2->ys = &yplot;
	series2->linearInterpolation = true;
	series2->lineType = toVector(L"solid");
	series2->lineThickness = 2;
	series2->color = CreateRGBColor(.1, .5, .8);


 
    // It is desired to find a curve (model function) of the form
    FunctionType function = [](const nc::NdArray<double>& coordinates, const nc::NdArray<double>& betas) -> double
    {
        const double s = coordinates.at(0);
        const double beta1 = betas.at(0);
        const double beta2 = betas.at(1);
 
        return (beta1 * s) / (beta2 + s);
    };
 
    // partial derivative of function with respect to beta1
    FunctionType delFdelBeta1 = [](const nc::NdArray<double>& coordinates, const nc::NdArray<double>& betas) -> double
    {
        const double s = coordinates.at(0);
        const double beta2 = betas.at(1);
 
        return s / (beta2 + s);
    };
 
    // partial derivative of function with respect to beta2
    FunctionType delFdelBeta2 = [](const nc::NdArray<double>& coordinates, const nc::NdArray<double>& betas) -> double
    {
        const double s = coordinates.at(0);
        const double beta1 = betas.at(0);
        const double beta2 = betas.at(1);
 
        return -(beta1 * s) / nc::square(beta2 + s);
    };
 
    // starting with the initial estimates of beta1Guess and beta2Guess and calculating after 5 iterations
    const nc::uint32 numIterations = 5;
    const double beta1Guess = 0.9;
    const double beta2Guess = 0.2;
 
#ifdef __cpp_structured_bindings
    auto [betas, rms] = nc::linalg::gaussNewtonNlls(numIterations, sMeasured.transpose(), rateMeasured,
        function, {delFdelBeta1, delFdelBeta2}, beta1Guess, beta2Guess);
#else
    auto results = nc::linalg::gaussNewtonNlls(numIterations, sMeasured.transpose(), rateMeasured,
        function, {delFdelBeta1, delFdelBeta2}, beta1Guess, beta2Guess);
    auto& betas = results.first;
    auto& rms = results.second;
    
#endif
 
    std::cout << "==========Wikipedia Example==========\n";
    std::cout << "beta values = " << betas;
    std::cout << "RMS = " << rms << '\n';
    
    std::vector<double> pred;
    for(auto i: x){
        pred.push_back(betas[0]*i/(betas[1]+i));
    }
    
    ScatterPlotSeries *series3 = GetDefaultScatterPlotSeriesSettings();
	series3->xs = &x;
	series3->ys = &pred;
	series3->linearInterpolation = true;
	series3->lineType = toVector(L"solid");
	series3->lineThickness = 2;
	series3->color = CreateRGBColor(.9, .5, .3);

	ScatterPlotSettings *settings = GetDefaultScatterPlotSettings();
	settings->width = 600;
	settings->height = 400;
	settings->autoBoundaries = true;
	settings->autoPadding = true;
	settings->title = toVector(L"mx+b Regression vs b1*x/(b2+x) Regression");
	settings->xLabel = toVector(L"substrate concentration");
	settings->yLabel = toVector(L"reaction rate");
	settings->scatterPlotSeries->push_back(series);
	settings->scatterPlotSeries->push_back(series2);
        settings->scatterPlotSeries->push_back(series3);

	DrawScatterPlotFromSettings(imageReference, settings, NULL);
	std::vector<double> *pngdata = ConvertToPNG(imageReference->image);
	WriteToFile(pngdata, "regs.png");
	DeleteImage(imageReference->image);
}
 
void exponentialExample()
{
    // United States population (in millions) and the corresponding year:
    nc::NdArray<double> year = nc::arange<double>(1.0, 9.0); // just use time points rather than the year
    nc::NdArray<double> population = { 8.3, 11.0, 14.7, 19.7, 26.7, 35.2, 44.4, 55.9 };
 
    // It is desired to find a curve (model function) of the form
    FunctionType exponentialFunction = [](const nc::NdArray<double>& coordinates, const nc::NdArray<double>& betas) -> double
    {
        const double t = coordinates.at(0);
        const double beta1 = betas.at(0);
        const double beta2 = betas.at(1);
 
        return beta1 * nc::exp(beta2 * t);
    };
 
    // partial derivative of function with respect to beta1
    FunctionType delFdelBeta1 = [](const nc::NdArray<double>& coordinates, const nc::NdArray<double>& betas) -> double
    {
        const double t = coordinates.at(0);
        const double beta2 = betas.at(1);
 
        return nc::exp(beta2 * t);
    };
 
    // partial derivative of function with respect to beta2
    FunctionType delFdelBeta2 = [](const nc::NdArray<double>& coordinates, const nc::NdArray<double>& betas) -> double
    {
        const double t = coordinates.at(0);
        const double beta1 = betas.at(0);
        const double beta2 = betas.at(1);
 
        return beta1 * t * nc::exp(beta2 * t);
    };
 
    // starting with the initial estimates of beta1Guess and beta2Guess and calculating after 5 iterations
    const nc::uint32 numIterations = 5;
    const double beta1Guess = 6.0;
    const double beta2Guess = 0.3;
 
#ifdef __cpp_structured_bindings
    auto [betas, rms] = nc::linalg::gaussNewtonNlls(numIterations, year.transpose(), population,
        exponentialFunction, {delFdelBeta1, delFdelBeta2}, beta1Guess, beta2Guess);
#else
    auto results = nc::linalg::gaussNewtonNlls(numIterations, year.transpose(), population,
        exponentialFunction, {delFdelBeta1, delFdelBeta2}, beta1Guess, beta2Guess);
    auto& betas = results.first;
    auto& rms = results.second;
#endif
 
    std::cout << "==========Exponential Population Example==========\n";
    std::cout << "beta values = " << betas;
    std::cout << "RMS = " << rms << '\n';
}
 
int main()
{   
    wikipediaExample();
    exponentialExample();
}
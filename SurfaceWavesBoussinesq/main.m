//
//  main.m
//  Boussinesq1DBathymetry
//
//  Created by Jeffrey J. Early on 2/19/13.
//
//

#import <Foundation/Foundation.h>
#import <GLNumericalModelingKit/GLNumericalModelingKit.h>

int main(int argc, const char * argv[])
{
    
    @autoreleasepool {
        
        /************************************************************************************************/
        /*		Typical user adjustable parameters															*/
        /************************************************************************************************/
        
        GLFloat significantWaveHeight = 1.0; // peak-to-peak amplitude of the wave, in meters
        GLFloat waveLength = 150; // wavelength of the sinusoid, in meters
        GLFloat sandbarHeight = 2; // sandbar height, in meters
        GLFloat waterDepthAtSandbar = 2; // depth of the water at the sandbar, in meters
        GLFloat runupToSlopeBreak = 1.0e3; // distance from the domain edge (where the forcing occurs) to the slope break
        GLFloat slopeBreakToSandbar = 5.0e3; // distance from the slope break, to the top of the sandbar
        GLFloat nPoints = 1024; // number of points used to discretize the domain.
        
        GLFloat totalTime = 500; // in seconds
        GLFloat outputInterval = 1; // in seconds
        NSURL *outputFile = [NSURL URLWithString: @"/Users/jearly/Desktop/BoussinesqWave.nc"];
		      
        /************************************************************************************************/
        /*		Compute various problem parameters														*/
        /************************************************************************************************/
        
        GLFloat g = 9.81;
        GLFloat beta = 0.537; // depth, relative to local bottom depth, where u will be solved. 0.537 corresponds to Nwogu's optimal value. 1-sqrt(1./3.);
        
        GLFloat shelfSlope = 5e-3; // Typical values range from 4-6 meters per kilometer
        GLFloat sandbarPosition = -500.0;
        GLFloat sandbarWidth = 100;
        
        GLFloat L_sponge = 4*waveLength/sqrt(2*M_PI); // This is the factor that goes in the decaying expontial for the sponge
        GLFloat spongeWidth = sqrt(4*L_sponge*L_sponge*log(10)); // This is where the sponge falls off to to 10^{-4} of its maximum, and where we will set 0 in the domain.
        
        GLFloat domainWidth = spongeWidth + fabs(sandbarPosition) + fabs(slopeBreakToSandbar) + fabs(runupToSlopeBreak);
        GLFloat domainMin = -(fabs(sandbarPosition) + fabs(slopeBreakToSandbar) + fabs(runupToSlopeBreak));
        GLFloat slopeBreakPosition = -(fabs(sandbarPosition) + fabs(slopeBreakToSandbar));
        
        // The function y=mx+b defining the linear slope.
        GLFloat m = -shelfSlope;
        GLFloat b = sandbarHeight + waterDepthAtSandbar - m*sandbarPosition;
        
        GLFloat depth = m*slopeBreakPosition + b; // the initial depth of the left most part of the domain.
        
        GLFloat alpha = beta*beta/2 - beta; // Nwogu's alpha parameter
        
        // The k,omega that we'll be forcing at for sin(omega*t - k*x)
        GLFloat k = 2*M_PI/waveLength;
        GLFloat omega = k*sqrt(g*depth)*sqrt(1-(alpha+1./3.)*depth*depth*k*k)/sqrt(1-alpha*depth*depth*k*k);
        
        // These are the coefficients for the forcing
        GLFloat N0 = significantWaveHeight/2.;
        GLFloat U0 = N0*(omega/(depth*k))/(1.-(alpha+1./3.)*depth*depth*k*k);
        GLFloat F_eta = U0*k*depth*(1.-(alpha+1./3.)*depth*depth*k*k);
        GLFloat F_u = g*k*N0;
        
        [GLBasisTransformOperation setShouldAntialias: NO];
        
        NSLog(@"Initial depth (D) is %f meters, and the wavelength (L) is %f meters, (D/L)^2 = %f", depth, waveLength, pow(depth/waveLength,2.0));
        NSLog(@"Initial wave height (N) to depth (D) ratio is N/D=%f", N0/depth);
        NSLog(@"Wave length set to %f meters, with oscillation time of %f seconds.", waveLength, 2*M_PI/omega);
        NSLog(@"u is being evaluated at %.3f h, and alpha=%.3f", beta, alpha);
        
        /************************************************************************************************/
        /*		Define the problem dimensions															*/
        /************************************************************************************************/
        
        GLEquation *equation = [[GLEquation alloc] init];
        
        GLDimension *xDim =[[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints: nPoints domainMin: domainMin length:domainWidth];
        xDim.name = @"x";
        xDim.units = @"meters";
        
        GLMutableDimension *tDim = [[GLMutableDimension alloc] initWithPoints: @[@(0.0)]];
        tDim.name = @"time";
        tDim.units = @"seconds";
        
        // This is the dimension for the cosine basis.
        GLDimension *kDim = [[GLDimension alloc] initAsDimension: xDim transformedToBasis: kGLCosineBasis strictlyPositive: YES];
        
        GLFunction *x = [GLFunction functionOfRealTypeFromDimension: xDim withDimensions: @[xDim] forEquation: equation];
        
        NSLog(@"Resolution set to %f meters.", xDim.sampleInterval);
        
        
        /************************************************************************************************/
        /*		Create the bathymetery                                                                  */
        /************************************************************************************************/
        
        GLFunction *bathymetry = [GLFunction functionOfRealTypeFromDimension: xDim withDimensions: @[xDim] forEquation: equation];
        
        bathymetry.name = @"bathymetry";
        for (NSUInteger i=0; i<bathymetry.nDataPoints; i++) {
            if (x.pointerValue[i] < slopeBreakPosition) {
                bathymetry.pointerValue[i] = depth;
            } else if (x.pointerValue[i] < 0.0) {
                bathymetry.pointerValue[i] = m*(x.pointerValue[i]) + b - sandbarHeight*exp(- pow(( x.pointerValue[i] - sandbarPosition)/sandbarWidth,2.0));
            } else {
                bathymetry.pointerValue[i] = b;
            }
            
        }
        
        /************************************************************************************************/
        /*		Create the initial conditions															*/
        /************************************************************************************************/
        
        GLFunction *initialEta = [GLFunction functionOfRealTypeWithDimensions: @[xDim] forEquation: equation];
        [initialEta zero];
        GLFunction *initialU = [initialEta scalarMultiply: sqrt(g/depth)];
        
        // Create the region of wave generation.
        GLFloat sigma = 3*xDim.sampleInterval;
        //sigma = waveLength/sqrt(2*M_PI);
        GLFunction *taper = [[[[[x plus: @(-xDim.domainMin)] times: @(1.0/sigma)] pow: 2.0] negate] exponentiate];
        taper.name = @"wave_generator";
        GLFunction *inverseTaper = [[taper negate] plus: @(1.0)];
        
        // Create the region of wave dissipation
        GLFloat spongeScale = -3.0*log(10.0)*sqrt(g*b)/L_sponge; // The fastest signals should decay to 10^{-3} their original values.
        GLFunction *sponge = [[[[[[x plus: @(-(xDim.domainMin + xDim.domainLength))] times: @(1.0/L_sponge)] pow: 2.0] negate] exponentiate] times: @(spongeScale)];
        sponge.name = @"sponge";
        
        /************************************************************************************************/
        /*		Create a file to output data															*/
        /************************************************************************************************/
        
        // Now we create a mutable variable in order to record the evolution of the Gaussian.
        GLNetCDFFile *netcdfFile = [[GLNetCDFFile alloc] initWithURL: outputFile forEquation: equation overwriteExisting: YES];
        GLMutableVariable *sshHistory = [initialEta variableByAddingDimension: tDim];
        sshHistory.name = @"SSH";
        sshHistory = [netcdfFile addVariable: sshHistory];
        
        [netcdfFile addVariable: bathymetry];
        [netcdfFile addVariable: taper];
        [netcdfFile addVariable: sponge];
        
        /************************************************************************************************/
        /*		Estimate the time step size																*/
        /************************************************************************************************/
        
        CGFloat cfl = 0.5;
        GLFloat U = sqrt(g*depth);
        GLFloat timeStep = cfl * xDim.sampleInterval / U;
        
        /************************************************************************************************/
        /*		Create and cache the differential operators we will need								*/
        /************************************************************************************************/
        
        // Create the operator L = (I + diag{alpha*h*h}*C^{-1}*d_xx*C - diag{beta*h*h_xx})
        GLLinearTransform *dct = [GLLinearTransform discreteTransformFromDimension: xDim toBasis: kGLCosineBasis forEquation:equation];
        GLLinearTransform *idct = [GLLinearTransform discreteTransformFromDimension: dct.toDimensions[0] toBasis: kGLDeltaBasis forEquation: equation];
        
        GLLinearTransform *cos_ddx = [GLLinearTransform differentiationMatrixFromDimension: dct.toDimensions[0] forEquation: equation];
        GLLinearTransform *sin_ddx = [GLLinearTransform differentiationMatrixFromDimension: cos_ddx.toDimensions[0] forEquation:equation];
        GLLinearTransform *ddxx = [sin_ddx times: cos_ddx];
        
        GLLinearTransform *B = [GLLinearTransform transformOfType: kGLRealDataFormat withFromDimensions: @[xDim] toDimensions: @[xDim] inFormat: @[@(kGLDenseMatrixFormat)] forEquation: equation];
        GLVariable *bathFactor = [[bathymetry times: bathymetry] scalarMultiply: alpha];
        [bathFactor solve];
        [B setVariableAlongDiagonal: bathFactor];
        
        GLLinearTransform *identity = [GLLinearTransform transformOfType: kGLRealDataFormat withFromDimensions: @[xDim] toDimensions: @[xDim] inFormat: @[@(kGLDenseMatrixFormat)] forEquation: equation];
        GLVariable *ones = [GLVariable variableOfRealTypeWithDimensions: @[xDim] forEquation:equation];
        ones = [ones setValue: 1.0 atIndices: @":"];
        [ones solve];
        [identity setVariableAlongDiagonal: ones];
        
        GLLinearTransform *B2 = [GLLinearTransform transformOfType: kGLRealDataFormat withFromDimensions: @[xDim] toDimensions: @[xDim] inFormat: @[@(kGLDenseMatrixFormat)] forEquation: equation];
        GLVariable *bathFactor2 = [[bathymetry times: [bathymetry xx]] scalarMultiply: -beta];
        [bathFactor2 solve];
        [B2 setVariableAlongDiagonal: bathFactor2];
        
        GLLinearTransform *operator = [[identity plus: B2] plus: [[[B times: idct] times: ddxx] times: dct]];
        //GLLinearTransform *operator = [identity plus: [[[B times: idct] times: ddxx] times: dct]];
        GLLinearTransform *inverseOperator = [operator inverse];
        NSLog(@"Inverting the differentiation matrix. This could take a while...");
        [inverseOperator solve];
        
        GLVariable *initialY = [operator transform: initialU];
        
        [equation setDefaultDifferentiationBasis: @[@(kGLCosineHalfShiftBasis)] forOrder: 1];
        
        diffOperators = [equation defaultDifferentialOperatorPoolForVariable: x];
        
        // Create our damping operator
        GLFloat nu = pow(xDim.sampleInterval/M_PI, 2.0)/timeStep;
        GLSpectralDifferentialOperator *svv = [diffOperators spectralVanishingViscosityFilter];
        //[diffOperators setDifferentialOperator: [[diffOperators harmonicOperatorOfOrder: 1] scalarMultiply: nu] forName: @"damp"];
        [diffOperators setDifferentialOperator: [[[diffOperators harmonicOperatorOfOrder: 1] scalarMultiply: nu] multiply: svv] forName: @"damp"];
        
        /************************************************************************************************/
        /*		Create an integrator: dy/dt=f															*/
        /************************************************************************************************/
        
        //GLRungeKuttaOperation *integrator = [GLAdaptiveRungeKuttaOperation rungeKutta23AdvanceY: @[initialY, initialEta] stepSize: timeStep fFromTY:^(GLVariable *time, NSArray *yNew) {
        GLRungeKuttaOperation *integrator = [GLRungeKuttaOperation rungeKutta4AdvanceY: @[initialY, initialEta] stepSize: timeStep fFromTY:^(GLVariable *time, NSArray *yNew) {
            
            GLVariable *u = [inverseOperator transform: yNew[0]];
            GLVariable *eta = yNew[1];
            
            GLVariable *fU = [[[[eta x] spatialDomain] scalarMultiply: -g] minus: [u times: [u x]]];
            
            GLVariable *bathymetryFactor1 = [bathymetry plus: [[[bathymetry times: bathymetry] times: [bathymetry xx]] scalarMultiply: -beta+0.5]];
            GLVariable *bathymetryFactor2 = [[[bathymetry times: bathymetry] times: bathymetry] scalarMultiply: alpha+1./3.];
            GLVariable *fEta = [[[[[[eta plus: bathymetryFactor1] times: u] plus: [bathymetryFactor2 times: [u xx]]] x] negate] spatialDomain];
            //GLVariable *fEta = [[[[[eta plus: bathymetry] times: u] x] negate] spatialDomain];
            
            // Apply the sponge
            fU = [fU plus: [yNew[0] times: sponge]];
            fEta = [fEta plus: [yNew[1] times: sponge]];
            
            // Create the wave generator
            GLVariable *easein = [[[[[time scalarMultiply: omega/(5*2.0*M_PI)] negate] exponentiate] negate] scalarAdd:1.0];
            GLVariable *arg = [[time scalarMultiply: omega] minus: [x scalarMultiply: k]];
            GLVariable *etaBoundary = [[arg cos] multiply: [easein scalarMultiply: F_eta]];
            GLVariable *uBoundary = [[arg cos] multiply: [easein scalarMultiply: F_u]];
            
            //			fU = [fU plus: [uBoundary times: taper]];
            //            fEta = [fEta plus: [etaBoundary times: taper]];
            
            fU = [[uBoundary times: taper] plus: [fU times: inverseTaper]];
            fEta = [[etaBoundary times: taper] plus: [fEta times: inverseTaper]];
            
            //			fEta = [[etaBoundary times: taper] plus: [etaBoundary times: inverseTaper]];
            //			fU = [[uBoundary times: taper] plus: [uBoundary times: inverseTaper]];
            
            //			fEta = [[etaBoundary times: taper] plus: [fEta times: inverseTaper]];
            //			fEta = [[etaBoundary times: taper] plus: [etaBoundary times: inverseTaper]];
            
            //			fU = [[uBoundary times: taper] plus: [fU times: inverseTaper]];
            //			fU = [[uBoundary times: taper] plus: [uBoundary times: inverseTaper]];
            
            
            
            //			GLVariable *etaBoundary = [[[time scalarMultiply: omega] cos] multiply: [easein scalarMultiply: F_eta]];
            //			GLVariable *uBoundary = [[[time scalarMultiply: omega] cos] multiply: [easein scalarMultiply: F_u]];
            //			fEta = [fEta setVariableValue: etaBoundary atIndices: @"0"];
            //			fU = [fU setVariableValue: uBoundary atIndices: @"0"];
            
            // Damp
            fU = [fU plus: [[u diff:@"damp"] spatialDomain]];
            fEta = [fEta plus: [[eta diff:@"damp"] spatialDomain]];
            
            return @[fU, fEta];
        }];
        
        /************************************************************************************************/
        /*		Now iterate! Stop every day to write out some data.										*/
        /************************************************************************************************/
        
        for (GLFloat time = outputInterval; time < totalTime; time += outputInterval)
        {
            @autoreleasepool {
                GLVariable *eta = [integrator stepForwardToTime: time][1];
                if (!eta)	break;
                
                NSLog(@"Logging time: %f s, last step size: %f, next step size: %f.", (integrator.currentTime), integrator.lastStepSize, integrator.stepSize);
                
                // We're using spectral code, so it's possible (and is in fact the case) that the variable is not in the spatial domain.
                [tDim addPoint: @(integrator.currentTime)];
                [sshHistory concatenateWithLowerDimensionalVariable: eta alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
            }
        }
        
        NSLog(@"Close the NetCDF file and wrap up");
        
        [equation waitUntilAllOperationsAreFinished];
        
        // The NetCDF file may still be writing data. We need to make sure it finishes before we exit.
        [netcdfFile waitUntilAllOperationsAreFinished];
        [netcdfFile close];
    }
    return 0;
}

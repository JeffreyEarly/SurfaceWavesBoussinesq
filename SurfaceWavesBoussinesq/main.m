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
        GLFloat waterDepthAtSandbar = 3.5; // depth of the water at the sandbar, in meters
        GLFloat runupToSlopeBreak = 0.5e3; // distance from the domain edge (where the forcing occurs) to the slope break
        GLFloat slopeBreakToSandbar = 4.0e3; // distance from the slope break, to the top of the sandbar
        GLFloat nPoints = 4096; // number of points used to discretize the domain.
		
		// Note: nPoints isn't scaling well. It requires the cfl to be 0.05, whereas 2048 requires cfl=0.25, 1024 works for less, etc.
		// Obviously I've scaled something incorrectly--damping?
		
        GLFloat totalTime = 1000; // in seconds
        GLFloat outputInterval = 1; // in seconds
        NSURL *outputFile = [NSURL URLWithString: @"/Users/jearly/Desktop/BoussinesqWaveHD.nc"];
		      
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
        
        GLDimension *xDim =[[GLDimension alloc] initDimensionWithGrid: kGLInteriorGrid nPoints: nPoints domainMin: domainMin length:domainWidth];
        xDim.name = @"x";
        xDim.units = @"meters";
		
		GLDimension *tDim = [[GLDimension alloc] initDimensionWithGrid: kGLEndpointGrid nPoints: 1 + round(totalTime/outputInterval)  domainMin:0 length:totalTime];
		tDim.name = @"time";
		tDim.units = @"seconds";
		
        GLFunction *x = [GLFunction functionOfRealTypeFromDimension: xDim withDimensions: @[xDim] forEquation: equation];
        
        NSLog(@"Resolution set to %f meters.", xDim.sampleInterval);
        
        
        /************************************************************************************************/
        /*		Create the bathymetery                                                                  */
        /************************************************************************************************/
        
        GLFunction *bathymetry = [GLFunction functionOfRealTypeFromDimension: xDim withDimensions: @[xDim] forEquation: equation];
		
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
		for (NSUInteger i=0; i<taper.nDataPoints; i++) {
			if (taper.pointerValue[i] < 1e-8) {
				taper.pointerValue[i] = 0.0;
			}
		}
        GLFunction *inverseTaper = [[taper negate] plus: @(1.0)];
        
        // Create the region of wave dissipation
        GLFloat spongeScale = -3.0*log(10.0)*sqrt(g*b)/L_sponge; // The fastest signals should decay to 10^{-3} their original values.
        GLFunction *sponge = [[[[[[x plus: @(-(xDim.domainMin + xDim.domainLength - spongeWidth/2.0))] times: @(2.0/L_sponge)] pow: 2.0] negate] exponentiate] times: @(spongeScale)];
		for (NSUInteger i=0; i<sponge.nDataPoints; i++) {
			if (x.pointerValue[i] > spongeWidth/2.0) {
				sponge.pointerValue[i] = spongeScale;
			}
		}
        
        /************************************************************************************************/
        /*		Estimate the time step size																*/
        /************************************************************************************************/
        
        CGFloat cfl = 0.1;
        GLFloat U = sqrt(g*depth);
        GLFloat timeStep = cfl * xDim.sampleInterval / U;
        
        /************************************************************************************************/
        /*		Create and cache the differential operators we will need								*/
        /************************************************************************************************/
        
        // Create the operator L = (I + diag{alpha*h*h}*C^{-1}*d_xx*C - diag{beta*h*h_xx})
        GLLinearTransform *dct = [GLLinearTransform discreteTransformFromDimension: xDim toBasis: kGLCosineBasis forEquation:equation];
        GLLinearTransform *idct = [GLLinearTransform discreteTransformFromDimension: dct.toDimensions[0] toBasis: kGLDeltaBasis forEquation: equation];
        GLLinearTransform *ddxx = [GLLinearTransform differentialOperatorWithDerivatives: 2 fromDimension: dct.toDimensions[0] forEquation: equation];
        
        GLLinearTransform *B = [GLLinearTransform transformOfType: kGLRealDataFormat withFromDimensions: @[xDim] toDimensions: @[xDim] inFormat: @[@(kGLDenseMatrixFormat)] forEquation: equation matrix: nil];
        GLFunction *bathFactor = [[bathymetry times: bathymetry] scalarMultiply: alpha];
        [bathFactor solve];
        [B setVariableAlongDiagonal: bathFactor];
        
        GLLinearTransform *identity = [GLLinearTransform transformOfType: kGLRealDataFormat withFromDimensions: @[xDim] toDimensions: @[xDim] inFormat: @[@(kGLDenseMatrixFormat)] forEquation: equation matrix:^( NSUInteger *row, NSUInteger *col ) {
            return (GLFloatComplex) (row[0]==col[0] ? 1.0 : 0.0);
        }];
        
        GLLinearTransform *B2 = [GLLinearTransform transformOfType: kGLRealDataFormat withFromDimensions: @[xDim] toDimensions: @[xDim] inFormat: @[@(kGLDenseMatrixFormat)] forEquation: equation matrix: nil];
        GLFunction *bathFactor2 = [[bathymetry times: [bathymetry xx]] scalarMultiply: -beta];
        [bathFactor2 solve];
        [B2 setVariableAlongDiagonal: bathFactor2];
        
        GLLinearTransform *operator = [[identity plus: B2] plus: [[[B times: idct] times: ddxx] times: dct]];
        //GLLinearTransform *operator = [identity plus: [[[B times: idct] times: ddxx] times: dct]];
        GLLinearTransform *inverseOperator = [operator inverse];
        NSLog(@"Inverting the differentiation matrix. This could take a while...");
        [inverseOperator solve];
		
        GLVariable *initialY = [operator transform: initialU];
		
		GLFloat nu = pow(xDim.sampleInterval/M_PI, 2.0)/timeStep;
		NSArray *spectralDimensions = dct.toDimensions;
		GLLinearTransform *laplacian = [GLLinearTransform harmonicOperatorOfOrder: 1 fromDimensions: spectralDimensions forEquation: equation];
		GLLinearTransform *svv = [GLLinearTransform spectralVanishingViscosityFilterWithDimensions: spectralDimensions scaledForAntialiasing: NO forEquation: equation];
		GLLinearTransform *harmonicDamp = [[laplacian times: @(nu)] multiply: svv];
        
        /************************************************************************************************/
        /*		Create an integrator: dy/dt=f															*/
        /************************************************************************************************/
		
		FfromTYVector flux = ^(GLScalar *time, NSArray *yNew) {
			
			GLFunction *u = [inverseOperator transform: yNew[0]];
			GLFunction *eta = yNew[1];
			
			GLFunction *fU = [[[[eta x] spatialDomain] scalarMultiply: -g] minus: [u times: [u x]]];
			//GLFunction *fU = [[[eta x] spatialDomain] scalarMultiply: -g];
			
			GLFunction *bathymetryFactor1 = [bathymetry plus: [[[bathymetry times: bathymetry] times: [bathymetry xx]] scalarMultiply: -beta+0.5]];
			GLFunction *bathymetryFactor2 = [[[bathymetry times: bathymetry] times: bathymetry] scalarMultiply: alpha+1./3.];
			GLFunction *fEta = [[[[[[eta plus: bathymetryFactor1] times: u] plus: [bathymetryFactor2 times: [u xx]]] x] negate] spatialDomain];
			//GLFunction *fEta = [[[[[eta plus: bathymetry] times: u] x] negate] spatialDomain];
			
			// Apply the sponge
			fU = [fU plus: [yNew[0] times: sponge]];
			fEta = [fEta plus: [yNew[1] times: sponge]];
			
			// Create the wave generator
			GLScalar *easein = [[[[[time times: @(omega/(5*2.0*M_PI))] negate] exponentiate] negate] plus: @(1.0)];
			GLFunction *arg = [[time times: @(omega)] minus: [x times: @(k)]];
			GLFunction *etaBoundary = [[arg cos] times: [easein times: @(F_eta)]];
			GLFunction *uBoundary = [[arg cos] times: [easein times: @(F_u)]];
			
			//			fU = [fU plus: [uBoundary times: taper]];
			//            fEta = [fEta plus: [etaBoundary times: taper]];
			
			fU = [[taper times: uBoundary] plus: [fU times: inverseTaper]];
			fEta = [[taper times: etaBoundary] plus: [fEta times: inverseTaper]];
			
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
			
			//			fEta = [fEta setValue: 0.0 atIndices: @"900:end"];
			//			fU = [fU setValue: 0.0 atIndices: @"900:end"];
			
			// Damp
			fU = [fU plus: [[[u transformToBasis: @[@(kGLCosineBasis)]] differentiateWithOperator: harmonicDamp] spatialDomain]];
			fEta = [fEta plus: [[[eta transformToBasis: @[@(kGLCosineBasis)]] differentiateWithOperator: harmonicDamp] spatialDomain]];
			
			return @[fU, fEta];
		};
		
		GLRungeKuttaOperation *integrator = [GLAdaptiveRungeKuttaOperation rungeKutta23AdvanceY: @[initialY, initialEta] stepSize: timeStep fFromTY:flux];
        //GLRungeKuttaOperation *integrator = [GLRungeKuttaOperation rungeKutta4AdvanceY: @[initialY, initialEta] stepSize: timeStep fFromTY:flux];
        
        /************************************************************************************************/
        /*		Now iterate! Stop every day to write out some data.										*/
        /************************************************************************************************/
		
		/************************************************************************************************/
		/*		Create a file to output data															*/
		/************************************************************************************************/
		
		// Now we create a mutable variable in order to record the evolution of the Gaussian.
		GLNetCDFFile *netcdfFile = [[GLNetCDFFile alloc] initWithURL: outputFile forEquation: equation overwriteExisting: YES];
		
		bathymetry = [bathymetry scaleVariableBy: 1.0 withUnits: @"m" dimensionsBy: 1.0 units: @"m"];
		bathymetry.name = @"bathymetry";
		
		taper = [taper scaleVariableBy: 1.0 withUnits: @"" dimensionsBy: 1.0 units: @"m"];
		taper.name = @"wave_generator";
		
		sponge = [sponge  scaleVariableBy: 1.0 withUnits: @"" dimensionsBy: 1.0 units: @"m"];
		sponge.name = @"sponge";
		
		[netcdfFile addVariable: bathymetry];
		[netcdfFile addVariable: taper];
		[netcdfFile addVariable: sponge];
		
		integrator.shouldDisplayProgress = YES;
		[integrator integrateAlongDimension: tDim withTimeScale: 1.0 file: netcdfFile output: ^(GLScalar *t, NSArray *y) {
			//NSLog(@"Logging day: %f, step size: %f.", (qg.T_QG*rkint.currentTime/86400), rkint.lastStepSize*qg.T_QG);
			
			NSArray *fluxVariables = flux(t, y);
			
			NSMutableDictionary *scaledVariables = [NSMutableDictionary dictionary];
			
			scaledVariables[@"SSH"] = [y[1] scaleVariableBy: 1.0 withUnits: @"m" dimensionsBy: 1.0 units: @"m"];
			scaledVariables[@"fU"] = [fluxVariables[0] scaleVariableBy: 1.0 withUnits: @"" dimensionsBy: 1.0 units: @"m"];
			scaledVariables[@"fEta"] = [fluxVariables[1] scaleVariableBy: 1.0 withUnits: @"" dimensionsBy: 1.0 units: @"m"];
			
			return scaledVariables;
		}];
        
        NSLog(@"Close the NetCDF file and wrap up");
        
        [equation waitUntilAllOperationsAreFinished];
        
        // The NetCDF file may still be writing data. We need to make sure it finishes before we exit.
        [netcdfFile waitUntilAllOperationsAreFinished];
        [netcdfFile close];
    }
    return 0;
}

#ifndef __OPTIMIZATIONFN_H
#define __OPTIMIZATIONFN_H

#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_block.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
extern gsl_rng * r;

/*
 *  optimizationfns.h
 *  
 *
 *  Created by Brian O'Meara on 4/7/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#define maxATTRACTIONPARAM         20
#define minATTRACTIONPARAM         0.001


class OptimizationFnMultiModel
{
	friend class BROWNIE;
public:
    BROWNIE brownie;
	OptimizationFnMultiModel(	gsl_matrix *Matrix0, gsl_matrix *Matrix1, gsl_matrix *Matrix2, gsl_matrix *Matrix3, gsl_matrix *Matrix4, gsl_matrix *Matrix5, gsl_matrix *Matrix6, gsl_matrix *Matrix7, gsl_matrix *Matrix8, gsl_matrix *Matrix9, gsl_vector *Vector1, gsl_vector *Vector2, int maxiterations, double stoppingprecision, int randomstarts, double stepsize, bool detailedoutput);  //VCV, observed values, tip variance
	~OptimizationFnMultiModel();
	int maxiterations;
	double stoppingprecision;
	int randomstarts;
	double stepsize;
	bool detailedoutput;
	
	gsl_vector * GeneralOptimization(int ChosenModel);

	double GetLikelihoodWithGivenTipVarianceOneRatePerState(const gsl_vector * variables);
	static double GetLikelihoodWithGivenTipVarianceOneRatePerState_gsl( const gsl_vector * variables, void *obj);
	gsl_vector * OptimizeRateWithGivenTipVarianceOneRatePerState();

	double GetLikelihoodWithGivenTipVarianceDiffRateOnChange(const gsl_vector * variables);
	static double GetLikelihoodWithGivenTipVarianceDiffRateOnChange_gsl( const gsl_vector * variables, void *obj);
	
	
	double GetLikelihoodOUSM(const gsl_vector * variables);
	static double GetLikelihoodOUSM_gsl( const gsl_vector * variables, void *obj);
	
	double GetLikelihoodOUSM_FixedAttractionRate(const gsl_vector * variables);
	static double GetLikelihoodOUSM_FixedAttractionRate_gsl( const gsl_vector * variables, void *obj);
	
	double GetLikelihoodOUSM_OnlyVariablesAttractionRate(const gsl_vector * variables);
	static double GetLikelihoodOUSM_OnlyVariablesAttractionRate_gsl( const gsl_vector * variables, void *obj);

	double GetLikelihoodOUSM_AnalyticMeans(const gsl_vector * variables);
	static double GetLikelihoodOUSM_AnalyticMeans_gsl( const gsl_vector * variables, void *obj);

	
//	gsl_vector * OptimizeRateWithGivenTipVarianceDiffRateOnChange();

//	double GetLikelihoodWithGivenTipVarianceOUFixedMeanFixedBrownVarAttract(const gsl_vector * variables);
//	static double GetLikelihoodWithGivenTipVarianceOUFixedMeanFixedBrownVarAttract_gsl( const gsl_vector * variables, void *obj);
//	gsl_vector * OptimizeRateWithGivenTipVarianceOUFixedMeanFixedBrownVarAttract();

//	double GetLikelihoodWithGivenTipVarianceOUFixedMeanFixedBrownVarAttract(const gsl_vector * variables);
//	static double GetLikelihoodWithGivenTipVarianceOUFixedMeanFixedBrownVarAttract_gsl( const gsl_vector * variables, void *obj);
//	gsl_vector * OptimizeRateWithGivenTipVarianceOUFixedMeanFixedBrownVarAttract();	
	
//	double GetLikelihoodWithGivenTipVarianceOUVarMeanFixedBrownVarAttract(const gsl_vector * variables);
//	static double GetLikelihoodWithGivenTipVarianceOUVarMeanFixedBrownVarAttract_gsl( const gsl_vector * variables, void *obj);
//	gsl_vector * OptimizeRateWithGivenTipVarianceOUVarMeanFixedBrownVarAttract();		
	
//	double GetLikelihoodWithGivenTipVarianceOUVarMeanVarBrownVarAttract(const gsl_vector * variables);
//	static double GetLikelihoodWithGivenTipVarianceOUVarMeanVarBrownVarAttract_gsl( const gsl_vector * variables, void *obj);
//	gsl_vector * OptimizeRateWithGivenTipVarianceOUVarMeanVarBrownVarAttract();		
	
	
private:
		gsl_matrix *Matrix0;
		gsl_matrix *Matrix1;
		gsl_matrix *Matrix2;
		gsl_matrix *Matrix3;
		gsl_matrix *Matrix4;
		gsl_matrix *Matrix5;
		gsl_matrix *Matrix6;
		gsl_matrix *Matrix7;
		gsl_matrix *Matrix8;
		gsl_matrix *Matrix9;
	gsl_vector *Vector1;
	gsl_vector *Vector2;
	gsl_vector *fixedparams;
};



class OptimizationFn
{
	friend class BROWNIE;
public:
    BROWNIE brownie;
	OptimizationFn(	gsl_matrix *Matrix1, gsl_vector *Vector1, gsl_vector *Vector2, int maxiterations, double stoppingprecision, int randomstarts, double stepsize, bool detailedoutput);  //VCV, observed values, tip variance
	~OptimizationFn();
	double GetLikelihoodWithGivenTipVariance(const gsl_vector * variables);
	static double GetLikelihoodWithGivenTipVariance_gsl( const gsl_vector * variables, void *obj);
	gsl_vector * OptimizeRateWithGivenTipVariance();

	double GetLikelihoodWithOptimizedTipVariance(const gsl_vector * variables);
	static double GetLikelihoodWithOptimizedTipVariance_gsl( const gsl_vector * variables, void *obj);

	double GetLikelihoodUnderACDC(const gsl_vector * variables);
	static double GetLikelihoodUnderACDC_gsl( const gsl_vector * variables, void *obj);

	double GetLikelihoodUnderOU1(const gsl_vector * variables);
	static double GetLikelihoodUnderOU1_gsl( const gsl_vector * variables, void *obj);

	double GetLikelihoodUnderDelta(const gsl_vector * variables);
	static double GetLikelihoodUnderDelta_gsl( const gsl_vector * variables, void *obj);

	double GetLikelihoodUnderLambda(const gsl_vector * variables);
	static double GetLikelihoodUnderLambda_gsl( const gsl_vector * variables, void *obj);
	
	gsl_vector * OptimizeRateWithOptimizedTipVariance();
        gsl_vector * GeneralOptimization(int ChosenModel);
	int maxiterations;
	double stoppingprecision;
	int randomstarts;
	double stepsize;
        bool detailedoutput;
	
private:
	gsl_matrix *Matrix1;
	gsl_vector *Vector1;
	gsl_vector *Vector2;
};

class LindyFn
{
	friend class BROWNIE;
public:
    BROWNIE brownie;
	TreesBlock* trees;
    TaxaBlock* taxa;
    AssumptionsBlock* assumptions;
    CharactersBlock* characters;	
	LindyFn( int maxiterations, double stoppingprecision, int randomstarts, double stepsize, bool detailedoutput,TreesBlock* trees, TaxaBlock* taxa, AssumptionsBlock* assumptions, CharactersBlock* characters); 
	int maxiterations;
	double stoppingprecision;
	int randomstarts;
	double stepsize;
	bool detailedoutput;	
	static double GetLikelihoodUnderLindy2_gsl( const gsl_vector * variables, void *obj) ;
		double GetLikelihoodUnderLindy2(const gsl_vector * variables);
		static double GetLikelihoodUnderLindy1_gsl( const gsl_vector * variables, void *obj) ;
		double GetLikelihoodUnderLindy1(const gsl_vector * variables);
		gsl_vector * GeneralOptimization(int ChosenModel);
};


//This class is for optimizing rate parameters on a given tree, with given model assignments on the branches. Up to 10 different models can be used
//class OptimizationFnDiscreteRate
//{
//	friend class BROWNIE;
//public:
//    BROWNIE brownie;
//	OptimizationFnDiscreteRate(	Tree *BrownieTree, gsl_matrix *Matrix0, gsl_matrix *Matrix1, gsl_matrix *Matrix2, gsl_matrix *Matrix3, gsl_matrix *Matrix4, gsl_matrix *Matrix5, gsl_matrix *Matrix6, gsl_matrix *Matrix7, gsl_matrix *Matrix8, gsl_matrix *Matrix9,  gsl_vector *Vector1, int maxiterations, double stoppingprecision, int randomstarts, double stepsize, bool detailedoutput, int statefrequencychoice);  
//	~OptimizationFnDiscreteRate();
//	int maxiterations;
//	double stoppingprecision;
//	int randomstarts;
//	double stepsize;
//	bool detailedoutput;
//	int statefrequencychoice;
	
//	gsl_vector * GeneralOptimization(int ChosenModel);
	
//	double GetLikelihoodWithGivenTipVarianceOneRatePerState(const gsl_vector * variables);
//	static double GetLikelihoodWithGivenTipVarianceOneRatePerState_gsl( const gsl_vector * variables, void *obj);
//	gsl_vector * OptimizeRateWithGivenTipVarianceOneRatePerState();
	
//	double GetLikelihoodWithGivenTipVarianceDiffRateOnChange(const gsl_vector * variables);
//	static double GetLikelihoodWithGivenTipVarianceDiffRateOnChange_gsl( const gsl_vector * variables, void *obj);
	
	
//	double GetLikelihoodOnTree(const gsl_vector * variables);
//	static double GetLikelihoodOnTree_gsl( const gsl_vector * variables, void *obj);
	

//private:
//		Tree *BrownieTree;
//	gsl_matrix *Matrix0;
//	gsl_matrix *Matrix1;
//	gsl_matrix *Matrix2;
//	gsl_matrix *Matrix3;
//	gsl_matrix *Matrix4;
//	gsl_matrix *Matrix5;
//	gsl_matrix *Matrix6;
//	gsl_matrix *Matrix7;
//	gsl_matrix *Matrix8;
//	gsl_matrix *Matrix9;
//	gsl_vector *Vector1;
//};


#endif
	
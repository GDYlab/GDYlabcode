This folder contains codes used for generating figures in the schema-based learning project:
These codes were implemented on Matlab 2020b and Ubuntu 18.04 platform

Fig1
~/Fig1/Fig1LearningCurves.m
Requires data from
~/Helpers/ChosenArmsTraining.m
~/Helpers/ChosenArmsExperiments.m
~/Helpers/GeneratePerformanceMatExperiments.m

Fig2
~/Fig2/Fig2Generalization.m
Calculate RSM matrices 
~/Helpers/CalcRSMSimilaritiesTemporalEachSessWithShufflesV2.m
~/Helpers/CalcExampleRSMPlotData.m
Needs Function
~/Helpers/cl2mat.m



Fig3
~/Fig3/Fig3Prediction.m
Requires function
~/Helpers/BayesDecode_FramesWithPlot.m
Requires data from
~/Helpers/CalcTimeSequencesInstantaneousError.m
~/Helpers/CalcStartArmTemporalDecodingNextTargetRTURT.m
~/Helpers/CalcStartArmDecodingTemporalEachTrialWithShuff.m
~/Helpers/CalcStartArmDecodingSpatialEachTrialWithShuff.m
~/Helpers/CalcStartArmSpatialDecodingNextTargetAllTrials50.m

Fig4
~/Fig4/Fig4SleepNetwork.m
~/Helpers/DecodeFramesConcatPriorBaselineAssimilated.m
~/Helpers/FindFramesSleep.m
~/Helpers/CalcLFPBestTet.m
~/Helpers/CalcReplayConcatPriorCombSleep.m
~/Helpers/CalcSleepRipples.m
~/Helpers/CalcTrackProbLearningEffectOnSleepFramesHPC.m
~/Helpers/CalcCorrMatForCuePlaceCoactivityInSleep.m
~/Helpers/CalcTimeSequenceConcatSleep.m
~/Helpers/CalcGoalProbLearningEffectsOnSleepFramesNewVsOld.m
~/Helpers/CalcTrackProbLearningEffectsOnSleepFramesNewVsOld.m
Needs Function
~/Helpers/CellAssembliesCorrMatSleepFrames.m


Fig5
~/Fig5/Fig5PFCCoordination.m
Requires Data from
~/Helpers/CellAssemblyFromAwakeRun.m
~/Helpers/CalcGoalFiringRatesCellAssemblies.m
~/Helpers/CalcGoalGenralizingPFCCAShuffleV2.m
~/Helpers/CalcRepresentationalMatrixForAwakeCellAssembliesGoalForSleepActivations.m
~/Helpers/CalcAwakeCellAssemblyPeriRippleActivationClassification.m
~/Helpers/CalcAwakeCAActivationInSleep25.m
~/Helpers/CalcMembershipCueCellsinPFCAssemblies.m
~/Helpers/CalcFiringRatesSleep.m
~/Helpers/CalcSingleCellPeriRippleActivation.m
~/Helpers/CalcFlavorCellsRippleActivityForFinalPlot.m
~/Helpers/CalcGoalProbConsolidationEffectNewVsOld.m
Needs function
~/Helpers/CellAssembliesToFiringRates.m

To get access to the formatted data, please contact the Dragoi's lab: george.dragoi@yale.edu

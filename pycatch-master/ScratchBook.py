import os
snapshot_number = str(1).zfill(2)
dir_name = os.path.join('./h' + snapshot_number + '/')
print(dir_name)

# for k, step in enumerate(np.arange(100, cfg.numberOfTimeSteps + 100, 100)):
#     snapshot_number = str(k).zfill(2)
#
#     if cfg.filtering:
#         if k > 0:
#             dir_name = os.path.join('./h' + snapshot_number + '/')
#             # if os.path.isdir(dir_name):
#             #     os.rmdir(dir_name)
#             os.rename('./1/', dir_name)
#
#         myModel = CatchmentModel()
#         dynamicModel = pcrfw.DynamicFramework(myModel, cfg.numberOfTimeSteps)
#         mcModel = pcrfw.MonteCarloFramework(dynamicModel, cfg.nrOfSamples)
#         mcModel.setForkSamples(True, 10)
#         # pfModel = SequentialImportanceResamplingFramework(mcModel)
#         pfModel = pcrfw.ResidualResamplingFramework(mcModel)
#         filterTimestepsNoSelection = range(3750, cfg.numberOfTimeSteps + 1, 25)
#         periodsToExclude = [
#             [2617, 2976],
#             [3649, 3689],
#             [4173, 4416],
#             [4046, 4366],
#             [5281, 6075]
#         ]
#         filterTimesteps = generalfunctions.removePeriodsFromAListOfTimesteps(filterTimestepsNoSelection,
#                                                                              periodsToExclude)
#         pfModel.setFilterTimesteps(filterTimesteps)
#         pfModel.run()
#
#     else:
#         if k > 0:
#             dir_name = os.path.join('./h' + snapshot_number + '/')
#             # if os.path.isdir(dir_name):
#             #     os.rmdir(dir_name)
#             os.rename('./1/', dir_name) # Rename raises an error if dir_name already exists,
#                                         # use lines above to remove dir_name if already present.
#
#         myModel = CatchmentModel()
#         dynamicModel = pcrfw.DynamicFramework(myModel, cfg.numberOfTimeSteps)
#         mcModel = pcrfw.MonteCarloFramework(dynamicModel, cfg.nrOfSamples)
#         mcModel.setForkSamples(True, 10)
#         mcModel.run()
#         dynamicModel.run()
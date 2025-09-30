import json

default_output = {
    "function": "output",
    "outputPath": "",
    "baseName": "fit",
    "printTrajectory": True,
    "verbose": True
}

default_optimizer_setting = {
    "function": "optimizerSetting",
    "lambdaL2": 1e-7,  # 0.000001,
    "pseudocount": 0,  # 20,
    "expBound": 40,
    "nThreads": 4,
    "nRetries": 3,
    "likelihoodThreshold": 0,  # 0.0002

}

default_lbfgsSettings = {
    "function": "lbfgsSettings",
    "memory": 100,
    "maxIters": 500,
    "convergence": 1e-7
}

default_setAlphabet = {
    "function": "setAlphabet",
    "letterOrder": "ACGT",
    "letterComplement": "C-G,A-T",
}

# no defaults in documentation, spec for rhoGamma models
default_enrichmentModelSeed = {
    "function": "enrichmentModelSeed",
    "rho": [1],
    "gamma": [0],
}

default_enrichmentModelConstraints = {
    "function": "enrichmentModelConstraints",
    "fitRho": False,
    "fitGamma": False,
    "roundSpecificRho": True,
    "roundSpecificGamma": True,
    "trySaturation": False
}


def generate_SMiLE_seq_configuration(
        count_table_path,
        variable_region_length, include_ns=True,
        gzipped="infer",
        right_flank="",
        left_flank="",
        binding_modes=2, binding_mode_size=12, binding_mode_flank=5, binding_mode_constraints=None
):
    """
    Generates a configuration for analysis of a classic SMiLE-seq experiment.
    :param include_ns: whether to include non-specific binding
    :param count_table_path: path to count table file
    :param variable_region_length: size of the variable sequence in the probe
    :param gzipped: indication if count table file is gzipped. Can infer if path has the ".gz" suffix.
    :param right_flank: right flank sequence
    :param left_flank: left flank sequence
    :param binding_modes: number of binding modes to model
    :param binding_mode_size: size of the binding site, can be an int (same value will be used for every binding mode)
            or a list of sizes for each binding mode
    :param binding_mode_flank: size of the analyzed flank (same value will be used for every binding mode)  or a list
            of sizes for each binding mode
    :param binding_mode_constraints: dictionary containing additional constraints for a binding mode (same value will
            be used for every binding mode)  or a list of dictionaries for each binding mode
    :return: ProBoundConfiguration object
    """
    # __generate_config_smile()
    config = __generate_config_smile(count_table_path, variable_region_length, gzipped=gzipped, include_ns=include_ns,
                                     right_flank=right_flank, left_flank=left_flank, binding_modes=binding_modes,
                                     binding_mode_size=binding_mode_size, binding_mode_flank=binding_mode_flank,
                                     binding_mode_constraints=binding_mode_constraints
                                     )
    return config


# with alternative workflow
def generate_meSMiLE_seq_configuration(
        count_table_path,
        variable_region_length,
        gzipped="infer",
        meCpG_encoding="mg",
        right_flank="",
        left_flank="", include_ns=True,
        binding_modes=2, binding_mode_size=12, binding_mode_flank=5, binding_mode_constraints=None
):
    """
        Generates a configuration for analysis of a meSMiLE-seq experiment.
        :param meCpG_encoding: how methylated CpG is indicated, must be two letters
        :param count_table_path: path to count table file
        :param variable_region_length: size of the variable sequence in the probe
        :param gzipped: indication if count table file is gzipped. Can infer if path has the ".gz" suffix.
        :param right_flank: right flank sequence
        :param left_flank: left flank sequence
        :param binding_modes: number of binding modes to model
        :param binding_mode_size: size of the binding site, can be an int (same value will be used for every binding mode)
                or a list of sizes for each binding mode
        :param binding_mode_flank: size of the analyzed flank (same value will be used for every binding mode)  or a list
                of sizes for each binding mode
        :param binding_mode_constraints: dictionary containing additional constraints for a binding mode (same value will
                be used for every binding mode)  or a list of dictionaries for each binding mode
        :return: ProBoundConfiguration object
        """
    # generate config, add alphabet
    config = __generate_config_smile(count_table_path, variable_region_length, gzipped=gzipped, include_ns=include_ns,
                                     right_flank=right_flank, left_flank=left_flank, binding_modes=binding_modes,
                                     binding_mode_size=binding_mode_size, binding_mode_flank=binding_mode_flank,
                                     binding_mode_constraints=binding_mode_constraints
                                     )
    config.alter_set_alphabet(letterOrder="ACGT"+meCpG_encoding,
                              letterComplement="C-G,A-T,"+f"{meCpG_encoding[0]}-{meCpG_encoding[1]}")
    return config


def __generate_config_smile(
        count_table_path,
        variable_region_length, gzipped="infer",
        right_flank="", left_flank="", include_ns=True,
        binding_modes=2,
        binding_mode_size=12,
        binding_mode_flank=5,
        binding_mode_constraints=None
):
    config = ProBoundConfiguration()
    config.alter_optimizer_setting(
        lambdaL2=0.000001,
        pseudocount=20,
        likelihoodThreshold=0.0002
    )
    config.add_table(count_table_path, 2, variable_region_length, gzipped_input=gzipped,
                     right_flank=right_flank, left_flank=left_flank)
    config.add_SELEX()
    if include_ns:
        config.add_nonspecific_binding()

    # sizes -- list OR integer
    # flanks -- list OR integer
    # constraints -- list OR dict

    # binding site size
    bms = [{} for i in range(binding_modes)]
    if isinstance(binding_mode_size, int):
        for item in bms:
            item["size"] = binding_mode_size
    elif isinstance(binding_mode_size, list):
        for item, val in zip(bms, binding_mode_size):
            item["size"] = val
    else:
        raise ValueError("Binding mode size can only be int or list")

    # flank length
    if isinstance(binding_mode_flank, int):
        for item in bms:
            item["flankLength"] = binding_mode_flank
    elif isinstance(binding_mode_flank, list):
        for item, val in zip(bms, binding_mode_flank):
            item["flankLength"] = val
    else:
        raise ValueError("Binding mode flank length can only be int or list")

    if binding_mode_constraints is None:
        for item in bms:
            config.add_binding_mode(size=item["size"], flankLength=item["flankLength"],
                                    binding_mode_constraints={
                                                  "maxFlankLength": -1,
                                                  "maxSize": 18,
                                                  "fittingStages": [
                                                      {"optimizeFlankLength": True},
                                                      {"optimizeMotifShiftHeuristic": True},
                                                      {"optimizeSizeHeuristic": True},
                                                  ]})
    elif isinstance(binding_mode_constraints, dict):
        for item in bms:
            config.add_binding_mode(size=item["size"], flankLength=item["flankLength"],
                                    binding_mode_constraints=binding_mode_constraints.copy())
    elif isinstance(binding_mode_constraints, list):
        for item, constr in zip(bms, binding_mode_constraints):
            config.add_binding_mode(size=item["size"], flankLength=item["flankLength"],
                                    binding_mode_constraints=constr)
    else:
        raise ValueError("Binding mode constraints can only be dict or list")

    return config


class ProBoundConfiguration:
    def __load_values(self, storage, fields, values):
        for field, value in zip(fields, values):
            if value is not None:
                storage[field] = value

    def __init__(self):
        self.__output = default_output
        self.__optimizerSettings = None
        self.__lbfgsSetting = None
        self.__setAlphabet = None
        self.__enrichmentModelSeed = None
        self.__enrichmentModelConstraints = None

        # these are private -- access them via methods
        self.__tables = []
        self.__binding_modes = []
        self.__interactions = []
        self.__enrichmentModels = []

        self.__binding_mode_index = -1
        self.__interactions_index = -1

        self.__gzip_oper = {
            "infer": lambda x: "tsv.gz" if x.endswith(".gz") else "tsv",
            "gz": lambda x: "tsv.gz",
            "tsv": lambda x: "tsv"
        }

    def get_all_setings(self):
        """
        Return all settings.
        :return: list of settings attributes
        """
        return [
            self.__output,
            self.__optimizerSettings,
            self.__lbfgsSetting,
            self.__setAlphabet,
            self.__enrichmentModelSeed,
            self.__enrichmentModelConstraints,
            self.__binding_modes,
            self.__interactions,
            self.__tables,
            self.__enrichmentModels
        ]

    def alter_output(self, output_path=None, base_name=None, print_trajectory=None, verbose=None):
        """
        Alter the default settings of output. If an attribute is set as None, it will not be altered (with regards to
        the current state of settings). Parameters match with the ProBound documentation.
        """
        if self.__output is None:
            self.__output = default_output.copy()  # copying -- default values will not be changed
        self.__load_values(self.__output,
                           ["outputPath", "baseName", "printTrajectory", "verbose"],
                           [output_path, base_name, print_trajectory, verbose])

    def alter_optimizer_setting(self, lambdaL2=None, pseudocount=None, expBound=None, nThreads=None, nRetries=None,
                                likelihoodThreshold=None):
        """
        Alter the default settings of optimizer. If an attribute is set as None, it will not be altered (with regards to
        the current state of settings). Parameters match with the ProBound documentation.
        """
        if self.__optimizerSettings is None:
            self.__optimizerSettings = default_optimizer_setting.copy()
        self.__load_values(self.__optimizerSettings,
                           ["lambdaL2", "pseudocount", "expBound", "nThreads", "nRetries", "likelihoodThreshold"],
                           [lambdaL2, pseudocount, expBound, nThreads, nRetries, likelihoodThreshold])

    def alter_lbfgs_setting(self, memory=None, maxIters=None, convergence=None):
        """
        Alter the default settings of lbfgs. If an attribute is set as None, it will not be altered (with regards to
        the current state of settings). Parameters match with the ProBound documentation.
        """
        if self.__lbfgsSetting is None:
            self.__lbfgsSetting = default_lbfgsSettings.copy()
        self.__load_values(self.__lbfgsSetting,
                           ["memory", "maxIters", "convergence"],
                           [memory, maxIters, convergence]
                           )

    def alter_set_alphabet(self, letterOrder=None, letterComplement=None):
        """
        Alter the default settings of the alphabet. If an attribute is set as None, it will not be altered (with
        regards to the current state of settings). Parameters match with the ProBound documentation.
        """
        if self.__setAlphabet is None:
            self.__setAlphabet = default_setAlphabet.copy()
        self.__load_values(self.__setAlphabet,
                           ["letterOrder", "letterComplement"],
                           [letterOrder, letterComplement]
                           )

    def alter_enrichment_model_seed(self, rho=None, gamma=None):
        """
        Alter the default settings of enrichment model seed (only rho gamma models). If an attribute is set as None,
        it will not be altered (with regards to the current state of settings). Parameters match with the ProBound
        documentation.
        """
        if self.__enrichmentModelSeed is None:
            self.__enrichmentModelSeed = default_enrichmentModelSeed.copy()
        self.__load_values(self.__enrichmentModelSeed,
                           ["rho", "gamma"],
                           [rho, gamma]
                           )

    def alter_enrichment_model_constraints(self, fit_rho=None, fitGamma=None,
                                           roundSpecificRho=None, roundSpecificGamma=None, trySaturation=None
                                           ):
        """
        Alter the default settings of enrichment model constraints(only rho gamma models). If an attribute is set as
        None, it will not be altered (with regards to the current state of settings). Parameters match with the
        ProBound documentation.
        """
        if self.__enrichmentModelConstraints is None:
            self.__enrichmentModelConstraints = default_enrichmentModelConstraints.copy()
        self.__load_values(self.__enrichmentModelConstraints,
                           ["fit_rho", "fitGamma", "roundSpecificRho", "roundSpecificGamma", "trySaturation"],
                           [fit_rho, fitGamma, roundSpecificRho, roundSpecificGamma, trySaturation]
                           )

    def add_table(self,
                  count_table_path,
                  nColumns,
                  variable_region_length,
                  gzipped_input="infer",
                  right_flank="",
                  left_flank="",
                  modeled_columns=None,
                  transliterate=None,
                  ):
        """
        Add a count table file to the config.
        :param count_table_path: path to the file
        :param nColumns: number of columns in the count table
        :param variable_region_length: probe size
        :param gzipped_input: is input gzipped
        :param right_flank: right flank sequence
        :param left_flank: left flank sequence
        :param modeled_columns: columns to model, default: all
        :param transliterate: change a letter to another
        """
        table_spec = {
                "function": "addTable",
                "leftFlank": left_flank,
                "rightFlank": right_flank,
                "variableRegionLength": variable_region_length,
                "countTableFile": count_table_path,
                "nColumns": nColumns,
                "inputFileType": self.__gzip_oper[gzipped_input](count_table_path)
            }
        if transliterate is not None:
            table_spec["transliterate"] = transliterate
        if modeled_columns is not None:
            table_spec["modeledColumns"] = modeled_columns

        self.__tables.append(table_spec)

    def add_SELEX(self,
                  modelType="SELEX",
                  bindingModes=None,
                  bindingModeInteractions=None,
                  cumulativeEnrichment=True,
                  concentration=1,
                  bindingSaturation=False
                  ):
        """
        Add an experiment. Parameters correspond to the ProBound arguments.
        """
        if bindingModes is None:
            bindingModes = [-1]
        if bindingModeInteractions is None:
            bindingModeInteractions = [-1]

        selex_spec = {
            "function": "addSELEX",
            "modelType": modelType,
            "bindingModes": bindingModes,
            "bindingModeInteractions": bindingModeInteractions,
            "cumulativeEnrichment": cumulativeEnrichment,
            "concentration": concentration,
            "bindingSaturation": bindingSaturation
        }
        self.__enrichmentModels.append(selex_spec)

    def add_nonspecific_binding(self):
        """
        Adds non-specific binding.
        """
        self.__binding_modes.append(
            {"function": "addNS"}
        )
        self.__binding_mode_index += 1

    def add_binding_mode(self, size,
                         flankLength=0, dinucleotideDistance=0, singleStrand=False, use_dinucleotides=False,
                         positionBias=False, binding_mode_constraints=None,  # for this binding mode
                         binding_mode_seed=None, symmetry=None):
        """
        Adds a specific binding mode. Parameters correspond to the ProBound arguments.
        The values of binding_mode_constraints, binding_mode_seed, symmetry can be set as dictionaries with parameters
        (as in the ProBound documentation) and there will be appropriated to the newly added binding mode.
        :param size: binding site size
        :param flankLength: size of the binding site flank
        :param dinucleotideDistance: maximum distance between the two letters of the dimer sequence features that are
            included in the model
        :param singleStrand: whether to only include forward strand
        :param use_dinucleotides: whether to use dinucleotides in motif inference. Sets dinucleotideDistance = 1
            (has precendence over explicit setting of dinucleotideDistance)
        :param positionBias: whether to include position bias
        :param binding_mode_constraints: the constraints on and the optimization process of binding mode parameters,
            a dictionary with appropriate keys (see ProBound manual)
        :param binding_mode_seed: the binding mode parameter seeds, a dictionary with appropriate keys
            (see ProBound manual)
        :param symmetry: symmetry status of a binding mode, a dictionary with appropriate keys (see ProBound manual)
        """
        if use_dinucleotides:
            dinucleotideDistance = 1
        self.__binding_modes.append(
              {
                "function": "addBindingMode",
                "size": size,
                "flankLength": flankLength,
                "dinucleotideDistance": dinucleotideDistance,
                "singleStrand": singleStrand,
                "positionBias": positionBias
              }
        )
        self.__binding_mode_index += 1
        # get the index of the just added binding mode
        bm_index = self.__binding_mode_index
        if binding_mode_constraints is not None:
            specs = {"function": "bindingModeConstraints", "index": bm_index}
            for k, val in binding_mode_constraints.items():
                specs[k] = val
            self.__binding_modes.append(specs)
        if binding_mode_seed is not None:
            specs = {"function": "bindingModeSeed", "index": bm_index}
            for k, val in binding_mode_seed.items():
                specs[k] = val
            self.__binding_modes.append(specs)
        if symmetry is not None:
            specs = {"function": "symmetry", "index": bm_index}
            for k, val in symmetry.items():
                specs[k] = val
            self.__binding_modes.append(specs)

    def get_binding_modes(self):
        """
        Get the list of the binding modes.
        :return:
        """
        return self.__binding_modes

    def add_interaction(self,
                        binding_mode1,
                        binding_mode2,
                        position_bias=False,
                        maxOverlap=0,
                        maxSpacing=-1,
                        interaction_constraints=None):
        """
        Add interaction between two binding modes.
        :param binding_mode1: index of interacting binding mode
        :param binding_mode2: index of interacting binding mode
        :param position_bias: see ProBound documentation
        :param maxOverlap: see ProBound documentation
        :param maxSpacing: see ProBound documentation
        :param interaction_constraints: dictionary of the interaction constraints to appropriate to this interaction
        """
        specs = {
            "function": "addInteraction",
            "bindingModes": [binding_mode1, binding_mode2],
            "positionBias": position_bias,
            "maxOverlap": maxOverlap,
            "maxSpacing": maxSpacing
        }
        self.__interactions.append(specs)
        interaction_index = self.__interactions_index
        if interaction_constraints is not None:
            specs = {"function": "interactionConstraints", "index": interaction_index}
            for k,v in interaction_constraints.items():
                specs[k] = v
            self.__interactions.append(specs)

    def get_interactions(self):
        """
        Get the list of interactions.
        """
        return self.__interactions

    def print_json(self, filename):
        """
        Write the config json.
        :param filename: path to output file
        :return:
        """
        first = True

        with open(filename, 'w') as json_writer:
            print("[", file=json_writer)
            for item in [self.__output, self.__optimizerSettings, self.__lbfgsSetting, self.__setAlphabet,
                         self.__enrichmentModelSeed, self.__enrichmentModelConstraints]:
                if item is not None:
                    if first:
                        first = False
                    else:
                        print(",", file=json_writer)
                    json.dump(item, json_writer)
            for col in [self.__tables, self.__enrichmentModels, self.__binding_modes, self.__interactions]:
                for item in col:
                    if first:
                        first = False
                    else:
                        print(",", file=json_writer)
                    json.dump(item, json_writer)
            print("\n]", file=json_writer)

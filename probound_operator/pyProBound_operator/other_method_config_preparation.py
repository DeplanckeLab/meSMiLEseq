from .input_config_preparation import ProBoundConfiguration

"""
Both of these are the most simplistic variants. This will be reworked once I have time for that.
"""

def generate_chip_config(count_table_path,
                         include_ns=True, probe_size=50, bm_size=12
                         ):
    config = ProBoundConfiguration()
    config.alter_optimizer_setting(
        lambdaL2=0.000001,
        pseudocount=20,
        likelihoodThreshold=0.0002
    )
    config.add_table(count_table_path,
                     2,
                     probe_size,
                     gzipped_input="infer",
                     left_flank="",
                     right_flank="")
    config.add_SELEX()
    if include_ns:
        config.add_nonspecific_binding()

    config.add_binding_mode(size=bm_size)
    # config.add_binding_mode(size=bm_size)
    return config


def generate_selex_confing(count_table_path, rounds, include_ns=True, probe_size=50, bm_size=12):
    config = ProBoundConfiguration()
    config.alter_optimizer_setting(
        lambdaL2=0.000001,
        pseudocount=20,
        likelihoodThreshold=0.0002
    )
    config.add_table(count_table_path,
                     rounds,
                     probe_size,
                     gzipped_input="infer",
                     left_flank="",
                     right_flank="")
    config.add_SELEX()
    if include_ns:
        config.add_nonspecific_binding()

    config.add_binding_mode(size=bm_size)
    # config.add_binding_mode(size=bm_size)
    return config

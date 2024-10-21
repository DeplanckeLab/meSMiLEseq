import json
import subprocess
import os

import logomaker
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd

jarpath = os.environ.get('PROBOUND_JAR_FULL_PATH')
print(f"Expecting ProBound .jar in {jarpath}")
probound_command = f"java -jar {jarpath}"

output_policies = {
    "rewrite": lambda config: delete_output(config),
    "append": lambda config: None  # does nothing
}


def delete_output(config_file):
    # read json filename
    with open(config_file) as json_handle:
        data = json.loads(json_handle.read())
    output_dict = [x for x in data if x["function"] == "output"][0]
    basename = output_dict["baseName"]
    outputname = f"{basename}.models.json"
    if os.path.exists(outputname):
        # this will ultimately fail if outputname is a directory but then you probably were doing
        # something really sketchy...
        os.remove(outputname)
        print(f"Existing file {outputname} was deleted.")


def run_probound(config_file,
                 full_config_file="tmp.fullconfig.json",
                 save_output="tmp.optimization.out",
                 cleanup_verbose=False, output_storage_policy="rewrite"):
    """
    Run the ProBound using Python.
    Be careful -- ProBound does not REWRITE the output file specified in the config (standard handling),
    it APPENDS to it. As a workaround, this func reads the config to infer the output file name and REMOVES IT BEFORE
    THE RUN!!!
    You can switch this behaviour off, but I STRONGLY ADVISE AGAINST IT.
    I honestly think this is a bug.
    :param output_storage_policy: what to do if an output model file with the same name already exists.
        Options: "rewrite", "append" (default of ProBound)
    :param config_file: path to the configuration file
    :param full_config_file: path to the full configuration file, will be deleted during the run
    :param save_output: path to file with ProBound command line messages
    :param cleanup_verbose: whether delete ProBound command line messages
    """
    output_storage_policy = output_storage_policy.lower()
    policy_implementation = output_policies[output_storage_policy]
    policy_implementation(config_file)

    my_env = os.environ.copy()
    my_env["PROBOUND_DIR"] = os.path.split(jarpath)[0]

    with open(full_config_file, mode="w") as fullconfig_handle:
        prep_process = subprocess.Popen(
            f"{probound_command} -b -c {config_file}",
            stdout=fullconfig_handle,
            stderr=subprocess.PIPE,
            text=True,
            shell=True,
            env=my_env
        )
    std_out, std_err = prep_process.communicate()
    if prep_process.returncode != 0:
        raise ValueError(f"Config file format error: {std_err}")

    with open(save_output, mode="w") as save_output_handle:
        run_process = subprocess.Popen(
            f"{probound_command} -c {full_config_file}",
            stdout=save_output_handle,
            stderr=subprocess.PIPE,
            text=True,
            shell=True,
            env=my_env
        )
    std_out, std_err = run_process.communicate()
    if run_process.returncode != 0:
        raise ValueError(f"Error in ProBound run: {std_err}")

    # cleanup
    os.remove(full_config_file)
    if cleanup_verbose:
        os.remove(save_output)


def __process_point(value_lines, sublevel_ident):
    result = value_lines
    particular_contents = []

    if len(value_lines) == 0:
        return [], {}

    i = 0
    line = value_lines[i]
    while line.split(" ")[0] != sublevel_ident:
        particular_contents.append(line)
        i += 1
        if i >= len(value_lines):
            break
        line = value_lines[i]

    sublevels = {}
    seen = None
    step_i = 0
    for line in value_lines[i:]:
        if line.split(" ")[0] == sublevel_ident:
            new_key = " ".join(line.split(" ")[1:])
            new_key = f"{step_i}: {new_key}"
            step_i += 1
            # print(line)
            if seen is not None:
                sublevels[key] = __process_point(seen, sublevel_ident + ">")
            key = new_key
            seen = []
        else:
            seen.append(line)
    if seen is not None:
        sublevels[key] = __process_point(seen, sublevel_ident + ">")

    return "\n".join(particular_contents), sublevels


def parse_probound_verbose(verbose_file):
    """
    A simple parser for command line messages. Probably not very useful.
    :param verbose_file: file with command line messages
    :return: dictionary with messages
    """
    with open(verbose_file) as probound_output:
        _, data = __process_point([line.strip("\n") for line in probound_output], ">")
    return data


def get_psam(output_json, dinucleotides=False, first_ns=True):
    """
    Basic data viewing of the calculated logos. Views the logos of the best (highest log likelihood) model.
    :param first_ns: whether first binding mode is non-specific and should be skipped (non-specific has no PSAM)
    :param dinucleotides: whether return dinucleotide psams as well
    :param output_json: path to the output json.
    :return: list of PSAM (pandas dataframe with letters/dimers as columns) for each binding mode
    if dinucleotides=True: list of PSAMs and list of dinucleotide PSAMs
    """
    model_data = []
    with open(output_json) as json_handle:  # they are liars this is not a json..........
        for line in json_handle:
            data = json.loads(line.strip("\n"))
            model_data.append(data)

    data = model_data[-1]

    alphabet = list(data["modelSettings"]["letterOrder"])
    coefs = [bm["mononucleotide"] for bm in data["coefficients"]["bindingModes"]]
    if first_ns:
        coefs = coefs[1:]

    # reorder coefs to tables w. r. to alphabet -- make pandas dfs
    psams = []
    for item in coefs:
        psam = np.reshape(item, (-1, len(alphabet)))
        psams.append(pd.DataFrame(psam, columns=alphabet))

    if dinucleotides:
        dipsams = []
        dicoefs = [bm["dinucleotide"] for bm in data["coefficients"]["bindingModes"]]
        dinucls = ["".join([x, y]) for x in alphabet for y in alphabet]

        for item in dicoefs:
            dipsam = np.reshape(item, (-1, len(dinucls)))
            dipsams.append(pd.DataFrame(dipsam, columns=dinucls))
        return psams, dipsams
    else:
        return psams


def view(psams):
    """
    Visualize every logo in the list of PSAMS. Very simple and only for convenience,
    use logomaker yourself if you want something fancier.
    :param psams: psams to plot
    :return: figure
    """
    fig, axes = plt.subplots(len(psams), ncols=1, figsize=(3 * len(psams), 6))
    if len(psams) == 1:
        logomaker.Logo(psams[0], ax=axes, color_scheme="classic")
    else:
        for psam, ax in zip(psams, axes):
            logomaker.Logo(psam, ax=ax, color_scheme="classic")

    plt.tight_layout()

    return fig, axes

import concurrent.futures
import copy
import multiprocessing
import os
import knime.extension as knext
import knime.types.chemistry as ktchem
import knime
import pandas as pd
from molvs import Standardizer
from rdkit import Chem, rdBase

root = os.path.dirname(os.path.abspath(__file__))

def is_mol(column):
    return column.ktype == knext.logical(Chem.Mol)

pharminfo = knext.category(
    path="/community",
    level_id="pharmacoinformatics",
    name="Pharmacoinformatics Research Group (UNIVIE)",
    description="Nodes created by the members of Pharmacoinformatics Research Group at the University of Vienna",
    icon="univie.png",
)
def process(mol, settings, debug=False):

    def is_nonorganic(fragment):
        return not any(a.GetAtomicNum() == 6 for a in fragment.GetAtoms())

    def contains_nonorg(fragment):
        return any(a.GetAtomicNum() not in [1, 6, 7, 8, 15, 16, 9, 17, 35, 53] for a in fragment.GetAtoms())

    s = Standardizer(max_tautomers=settings['max_tau'])
    functions = {
        "disconnect_metals": s.disconnect_metals,
        "remove_fragments": s.remove_fragments,
        "uncharge": s.uncharge,
        "normalize": s.normalize,
        "reionize": s.reionize,
        "stereo": {
            'KEEP': lambda x: x,
            'REMOVE': s.stereo_parent,
            'CLEAN': lambda x: (Chem.AssignStereochemistry(x, force=True, cleanIt=True), x)[1]
        },
        "largest_fragment": s.largest_fragment,
        "canonicalize_tautomer": s.canonicalize_tautomer,
        "is_inorganic": is_nonorganic,
        "contains_nonorganic": contains_nonorg,
    }

    def apply_with_info(m, f, e):
        if not m: return None, e
        r, err = None, ''
        try: r = f(m)
        except Exception as exc:
            err = str(exc)
        return (r, e + [err]) if err else (r, e)

    def apply(m, f):
        if not m:
           return None
        try:
            return  f(m)
        except:
            return None

    errors = []
    is_nonorganic_flag = contains_nonorg_flag = None

    for k, f in functions.items():
        if callable(f) and settings[k] and k not in ['is_inorganic', 'contains_nonorganic']:
            if debug:
                mol, errors = apply_with_info(mol, f, errors)
            else:
                mol = apply(mol, f)
        elif isinstance(f, dict):
            if debug:
                mol, errors = apply_with_info(mol, f[settings['stereo']], errors)
            else:
                mol = apply(mol, f[settings['stereo']])
        elif k == 'is_inorganic':
            is_nonorganic_flag = f(mol)
        elif k == 'contains_nonorganic':
            contains_nonorg_flag = f(mol)
        if not mol: break

    return (mol, is_nonorganic_flag, contains_nonorg_flag, "\n".join(errors)) if debug else (mol, is_nonorganic_flag, contains_nonorg_flag)


@knext.parameter_group("Univie Standardizer settings")
class Settings:
    molecule_column = knext.ColumnParameter(label="Molecules",description="The RDKit molecules column to proceed",column_filter=is_mol, port_index=0,include_none_column=False)
    class StereoOptions(knext.EnumParameterOptions):
        KEEP = ("No Changes", "Do nothing with stereochemistry.")
        REMOVE = ("Remove", "Remove stereochemistry at all.")
        CLEAN = ("Clean", "Clean and fix stereochemistry.")

    stereo = knext.EnumParameter(
        "Stereochemistry processing",
        "How to process stereochemistry.",
        StereoOptions.KEEP.name,
        StereoOptions,
    )

    debug = knext.BoolParameter("Debug", "Debug with information from RDkit core (slow, single-core)", False)
    sanitize = knext.BoolParameter("Sanitize", "Sanitize molecules (RDkit)", True)
    contains_nonorganic = knext.BoolParameter("Check Non-organic", "Check if the molecule contains non-organic fragments, adding additional field is_nonorganic to the final table", True)
    is_inorganic = knext.BoolParameter("Is Inorganic", "Check if the molecule is fully inorganic", True)
    remove_fragments = knext.BoolParameter("Remove Fragments", "Remove Fragments", True)
    uncharge = knext.BoolParameter("Uncharge", "Uncharge the molecules", True)
    disconnect_metals = knext.BoolParameter("Disconnect Metals", "Disconnect metal atoms", True)
    normalize = knext.BoolParameter("Normalize", "Normalize molecules", True)
    reionize = knext.BoolParameter("Reionize", "Reionize molecules", True)
    largest_fragment = knext.BoolParameter("Keep Largest Fragment", "Keep only the largest fragment", True)
    # enumerate_tautomers = knext.BoolParameter("Enumerate Tautomers", "Enumerate all possible tautomers", True)
    canonicalize_tautomer = knext.BoolParameter("Canonicalize Tautomer", "Convert tautomers to their canonical form",
                                                True)
    max_tau = knext.IntParameter("Max tautamers", "Maximum number of tautomer generated from a single molecule", 1000)

@knext.node(name="Standardizer", node_type=knext.NodeType.MANIPULATOR, icon_path="univiestand.ico", category=pharminfo)
@knext.input_table(name="Molecular Table", description="A dataset with molecular structures")
@knext.output_table(name="Standardized Molecules", description="The same datatable, plus standardized columns")
@knext.output_table(name="Errors", description="Table with molecules that caused errors")
class UniVieKnimeNode(knext.PythonNode):
    settings = Settings()

    def configure(self, config_context: knext.ConfigurationContext, input_schema: knext.Schema) -> knext.Schema:
        input_schema = input_schema.append(knext.Column(Chem.rdchem.Mol, "Standardized Molecules (UNIVIE)"))
        if self.settings.contains_nonorganic:
            input_schema = input_schema.append(knext.Column(knime.api.schema.bool_(),"Contains Inorganic"))
        if self.settings.is_inorganic:
            input_schema = input_schema.append(knext.Column(knime.api.schema.bool_(),"Is inorganic"))
        error_schema = knext.Schema([Chem.rdchem.Mol], names=["Erroneous Molecules"])
        if self.settings.debug:
            # input_schema = input_schema.append(knext.Column(knime.api.schema.string(),"RDkit Debug"))
            error_schema = error_schema.append(knext.Column(knime.api.schema.string(),"RDkit Error"))
        return input_schema, error_schema



    def execute(self, exec_context, input_1):
        df = input_1.to_pandas()
        if self.settings.molecule_column is None:
            raise ValueError("Please select a column with molecules!")

        settings = copy.copy(self.settings.__dict__["__parameters__"])
        errors = {}
        if self.settings.debug:
            debug_string_normal = {}
            debug_string_failed = {}
        results = {}
        is_inorganic = {}
        contains_nonorg = {}

        if not self.settings.debug:
            n_cores = multiprocessing.cpu_count()
            with concurrent.futures.ProcessPoolExecutor(max_workers=n_cores) as executor:
                futures = {executor.submit(process, mol,settings,False): index for index, mol in df[self.settings.molecule_column].items()}
                for future in concurrent.futures.as_completed(futures):
                    df_index = futures[future]
                    try:
                        results_row,is_inorganic_row,contains_nonorg_row = future.result()
                        if results_row is None:
                            errors[df_index] = df[self.settings.molecule_column][df_index]
                        else:
                            results[df_index] = results_row
                            is_inorganic[df_index] = is_inorganic_row
                            contains_nonorg[df_index] = contains_nonorg_row
                    except Exception as exc:
                        errors[df_index]  = df[self.settings.molecule_column][df_index]
        else:
            # If debug mode is enabled, process each molecule individually to capture error messages
            for index, mol in df[self.settings.molecule_column].items():
                try:

                    result, is_inorganic_flag, contains_nonorg_flag, errors_string = process(mol, settings, True)
                    if result is None:
                        errors[index] = df[self.settings.molecule_column][index]
                        debug_string_failed[index] = errors_string
                    else:
                        results[index] = result
                        is_inorganic[index] = is_inorganic_flag
                        contains_nonorg[index] = contains_nonorg_flag
                        debug_string_normal[index] = errors_string
                except Exception as exc:
                    # Save the error message and the molecule that caused it
                    errors[index] = df[self.settings.molecule_column][index]
                    debug_string_failed[index] = str("GLOBAL: " +exc.__class__.__name__ + ": " + str(exc))

        results_df = pd.DataFrame.from_dict(results, orient='index', columns=['Standardized Molecules (UNIVIE)'])
        if self.settings.contains_nonorganic:
            contains_nonorg_df = pd.DataFrame.from_dict(contains_nonorg, orient='index', columns=['Contains Inorganic'])
            results_df = results_df.merge(contains_nonorg_df, left_index=True, right_index=True)
        if self.settings.is_inorganic:
            is_inorganic_df = pd.DataFrame.from_dict(is_inorganic, orient='index', columns=['Is inorganic'])
            results_df = results_df.merge(is_inorganic_df, left_index=True, right_index=True)
        df = df.merge(results_df, left_index=True, right_index=True)
        errors_df = pd.DataFrame.from_dict(errors, orient='index', columns=['Erroneous Molecules'])
        if self.settings.debug:
            debug_string_failed_df = pd.DataFrame.from_dict(debug_string_failed, orient='index', columns=['RDkit Error'])
            errors_df = errors_df.merge(debug_string_failed_df, left_index=True, right_index=True)
        return knext.Table.from_pandas(df), knext.Table.from_pandas(errors_df)
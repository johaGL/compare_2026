import sys
import os
import logging
import click
import warnings


@click.command()
@click.argument('input_file', type=click.Path(exists=True))
@click.argument('data_dir', type=click.Path(file_okay=False, dir_okay=True))



def main(input_file, data_dir):
    """
    \b
    Obtain the metabolite resolved matrix only, per set type (WIKI and KEGG)
      Requires:
        - the `smpath` virtual environment
        - the scripts1/ already applied, that copies from '../../data/' to ../../compar_out_2026 all the required files (clean annotations, .h5ad and

    usage:
    cd scripts1
    python run.py ../the_data_list.tsv ../../compar_out_2026

    \b
    Positional arguments:

      INPUT_FILE: Relative or absolute path to the table in tab delimited format,
       where each dataset and user-defined parameters for its analysis are stored.

      DATA_DIR: Relative or absolute path, to the folder where the datasets are stored, e.g. data/.


    """
    click.echo(f'Input file: {input_file}')
    click.echo(f'Data directory: {data_dir}')
    # start
    try:
        if (input_file is not None) and (
            data_dir is not None
        ) :
            warnings.simplefilter('ignore', category=FutureWarning)
            from src.dispatch_tissues_slices import process_tissues
            process_tissues(input_file, data_dir)

    except Exception as e:
        logging.error(e)


if __name__ == '__main__':
    main()

set -e
# Clear the output directory except for the gene_list.txt
find batch_cloning_output -type f ! -name 'gene_list.txt' -delete
find batch_cloning_output -type d -empty -delete

python -m opencloning.batch_cloning.pombe.pombe_get_primers --genes batch_cloning_output/gene_list.txt
python -m opencloning.batch_cloning.pombe.pombe_clone --genes batch_cloning_output/gene_list.txt
python -m opencloning.batch_cloning.pombe.pombe_summary
python -m opencloning.batch_cloning.pombe.pombe_gather

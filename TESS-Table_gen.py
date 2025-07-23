# Function to format the LaTeX table
def format_latex_table(df):
    # Only keep necessary columns and reset index
    df = df[["KIC", "TIC", "R2FFT"]].reset_index(drop=True)

    # Split the dataframe into two halves
    half = (len(df) + 1) // 2
    first_half = df.iloc[:half].reset_index(drop=True)
    second_half = df.iloc[half:].reset_index(drop=True)

    # Pad second half with empty rows if needed
    while len(second_half) < len(first_half):
        second_half = second_half._append(pd.Series([None]*3, index=second_half.columns), ignore_index=True)

    # Create the LaTeX table
    table = r"""\begin{table}[h]
\centering
\renewcommand{\arraystretch}{1.5}
\begin{tabular}{|c|c|c||c|c|c|}
\hline
\textbf{KIC} & \textbf{TIC} & \textbf{$R^2_{FFT}$} & \textbf{KIC} & \textbf{TIC} & \textbf{$R^2_{FFT}$} \\
\hline
"""

    for i in range(len(first_half)):
        row1 = first_half.iloc[i]
        row2 = second_half.iloc[i]

        def fmt(val):
            if pd.isna(val):
                return ""
            if isinstance(val, float):
                return f"{val:.6f}"
            return str(val)

        table += f"{fmt(row1['KIC'])} & {fmt(row1['TIC'])} & {fmt(row1['R2FFT'])} & " \
                 f"{fmt(row2['KIC'])} & {fmt(row2['TIC'])} & {fmt(row2['R2FFT'])} \\\\\n\\hline\n"

    table += r"\end{tabular}" + "\n\\caption{KIC, TIC, and $R^2_{FFT}$ values.}\n\\end{table}"
    return table

import pandas as pd

# Load the CSV file
file_path = r"C:\Users\ahmed\Downloads\FINAL_DATA_TESS_COMP.csv"
df = pd.read_csv(file_path)

# Display the first few rows to inspect the structure
df.head()
latex_output = format_latex_table(df)
#latex_output[:1500]  # Preview the first part of the LaTeX output
# Save the LaTeX output to a .tex file
output_path = "kic_tic_r2fft_table.tex"
with open(output_path, "w") as f:
    f.write(latex_output)

#output_path

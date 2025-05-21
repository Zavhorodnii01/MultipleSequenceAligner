from multiple_sequence_aligner import MultipleSequenceAligner
import tkinter as tk
from tkinter import ttk
from PIL import Image, ImageDraw, ImageFont
import json
import os
from tkinter import filedialog, messagebox
import re

root = tk.Tk()
root.title("MSA")
root.geometry("1000x700")
root.configure(bg="#f4f4f4")

# Introductory Text
intro_label = tk.Label(root, text="This is an interactive example of"
    " Center Start Method used for Multiple Sequence Alignment",
    font=("Arial", 12), justify="center", wraplength=900)
intro_label.pack(pady=2)

# Style configuration
style = ttk.Style()
style.theme_use("default")
style.configure("TNotebook.Tab", font=('Helvetica', 12, 'bold'), padding=[20, 10], foreground="#555")
style.map("TNotebook.Tab",
              background=[("selected", "#ffffff")],
              foreground=[("selected", "#000")])
style.layout("TNotebook.Tab", [
        ('Notebook.tab', {
            'sticky': 'nswe',
            'children': [
                ('Notebook.padding', {
                    'side': 'top',
                    'sticky': 'nswe',
                    'children': [
                        ('Notebook.label', {'side': 'top', 'sticky': ''})
                    ]
                })
            ]
        })])
style.configure("TNotebook", tabposition='n')

header = tk.Label(root, text="MSA", font=("Helvetica", 22, "bold"), bg="#f4f4f4", fg="#4e8074")
header.pack(anchor="w", padx=10, pady=(10, 0))

notebook = ttk.Notebook(root)
notebook.pack(expand=True, fill="both", padx=10, pady=10)

# Initializes tabs names
tab_names = ["Input", "Parameters", "CLUSTAL Alignment", "FASTA Alignment", "AlignmentViewer"]
tabs = {}
for name in tab_names:
    frame = ttk.Frame(notebook, padding=20)
    notebook.add(frame, text=name)
    tabs[name] = frame

notebook.tab(2, state='disabled')
notebook.tab(3, state='disabled')
notebook.tab(4, state='disabled')

# Function for mouse scrolling
def bind_mousewheel(widget, text_widget):
    def _on_mousewheel(event):
        if event.state & 0x0001:
            text_widget.xview_scroll(-1 * int(event.delta / 120), "units")
        else:
            text_widget.yview_scroll(-1 * int(event.delta / 120), "units")

    widget.bind("<MouseWheel>", _on_mousewheel)
    widget.bind("<Shift-MouseWheel>", _on_mousewheel)
    widget.bind("<Button-4>", lambda e: text_widget.yview_scroll(-1, "units"))  # Linux scroll up
    widget.bind("<Button-5>", lambda e: text_widget.yview_scroll(1, "units"))   # Linux scroll down

# Sets the prompt
placeholder = "Enter protein/nucleotide sequences in FASTA format"
def setup_placeholder(text_field ):
    def on_key_press(event):
        current_text = text_field.get("1.0", "end-1c")
        if current_text == placeholder:
            text_field.delete("1.0", "end")
            text_field.config(fg="black")

    def on_focus_out(event):
        if text_field.get("1.0", "end-1c").strip() == "":
            text_field.insert("1.0", placeholder)
            text_field.config(fg="grey")

    text_field.insert("1.0", placeholder)
    text_field.config(fg="grey")
    text_field.bind("<Key>", on_key_press)
    text_field.bind("<FocusOut>", on_focus_out)

# loaded files
fasta_contents = []
def on_submit_btn_click(input_frame_dict, parameters_frame_dict):
    """
        Handles the submit button click event.

        - Updates the button label to "Resubmit"
        - Extracts and validates user input and parameters
        - Runs the multiple sequence alignment
        - Displays results in various formats
        - Switches to the results tab if successful

        Args:
            input_frame_dict (dict): Dictionary containing UI elements for input.
            parameters_frame_dict (dict): Dictionary containing UI elements for alignment parameters.
        """
    input_frame_dict['submit_btn'].config(text="Resubmit")

    global fasta_contents
    user_input = input_frame_dict['text_field'].get("1.0", "end-1c")

    try:
        all_user_input = extract_names_and_sequences("\n".join(fasta_contents) + user_input)
    except Exception as e:
        return

    if all_user_input:
        pass
    else:
        messagebox.showerror("Incorrect input","Incorrect Input, try again")

    def is_valid_number(s):
        return (len(s) > 0 and s.lstrip('-').isdigit() and
                (s.count('-') <= 1 and (s.startswith('-') or '-' not in s)))

    if (is_valid_number(parameters_frame_dict['match_score'].get()) == False
            or is_valid_number(parameters_frame_dict['mismatch_score'].get()) == False
            or is_valid_number(parameters_frame_dict['gap_score'].get()) == False
            or is_valid_number(parameters_frame_dict['match'].get()) == False
            or is_valid_number(parameters_frame_dict['substitution'].get()) == False
            or is_valid_number(parameters_frame_dict['gap'].get()) == False):
        messagebox.showerror("Invalid parameters", "Please enter valid numbers for match, mismatch and gap scores.")
        return

    match_score = int(parameters_frame_dict['match_score'].get())
    mismatch_score = int(parameters_frame_dict['mismatch_score'].get())
    gap_score = int(parameters_frame_dict['gap_score'].get())
    match = int(parameters_frame_dict['match'].get())
    substitution = int(parameters_frame_dict['substitution'].get())
    gap = int(parameters_frame_dict['gap'].get())
    aligned_sequences, score, statistics = (get_aligned_sequences_score_statistics
                                            (all_user_input, match_score, mismatch_score, gap_score, match,
                                             substitution, gap))
    print_result_in_clustal_format(aligned_sequences, score, statistics, parameters_frame_dict)
    print_result_in_fasta_format(aligned_sequences, score, statistics, parameters_frame_dict)
    print_result_in_alignment_viewer(aligned_sequences, score, statistics)

    notebook.tab(2, state='normal')
    notebook.tab(3, state='normal')
    notebook.tab(4, state='normal')

    # moves to the result tab
    notebook.select(2)

def on_reset_btn_click(input_frame_dict, parameters_frame_dict):
    """
       Handles the reset button click event.

       - Clears the input text field and file message
       - Resets the submit button text and input placeholders
       - Disables result tabs
       - Resets global FASTA content list
       - Resets all scoring parameters to default values

       Args:
           input_frame_dict (dict): Dictionary with UI elements for user input.
           parameters_frame_dict (dict): Dictionary with UI elements for alignment parameters.
       """

    input_frame_dict['text_field'].delete("1.0", "end")
    input_frame_dict['load_file_message'].config(text="")


    setup_placeholder(input_frame_dict['text_field'])
    global fasta_contents
    fasta_contents = []

    input_frame_dict['submit_btn'].config(text="Submit")

    notebook.tab(2, state='disabled')
    notebook.tab(3, state='disabled')
    notebook.tab(4, state='disabled')
    parameters_frame_dict['match_score'].delete(0, tk.END)
    parameters_frame_dict['match_score'].insert(0, "1")
    parameters_frame_dict['mismatch_score'].delete(0, tk.END)
    parameters_frame_dict['mismatch_score'].insert(0, "-1")
    parameters_frame_dict['gap_score'].delete(0, tk.END)
    parameters_frame_dict['gap_score'].insert(0, "-2")
    parameters_frame_dict['match'].delete(0, tk.END)
    parameters_frame_dict['match'].insert(0, "1")
    parameters_frame_dict['substitution'].delete(0, tk.END)
    parameters_frame_dict['substitution'].insert(0, "-1")
    parameters_frame_dict['gap'].delete(0, tk.END)
    parameters_frame_dict['gap'].insert(0, "-2")


def get_aligned_sequences_score_statistics(user_input, match_score, mismatch_score, gap_score, match, substitution, gap):
    """
        Performs multiple sequence alignment and returns the results.

        Args:
            user_input (list): List of (name, sequence) tuples.
            match_score (int): Score for matching characters.
            mismatch_score (int): Penalty for mismatched characters.
            gap_score (int): Penalty for introducing a gap.
            match (int): Weight for match in statistics.
            substitution (int): Weight for substitution in statistics.
            gap (int): Weight for gap in statistics.

        Returns:
            tuple: (aligned_sequences, alignment_score, statistics)
        """
    msa = MultipleSequenceAligner(user_input, match_score, mismatch_score, gap_score, match, substitution, gap)
    return msa.get_final_alignments(), msa.get_score(), msa.get_statistics()

def print_result_in_clustal_format(final_alignments, score, statistics, parameters_frame_dict):
    clustal_alignment = tabs["CLUSTAL Alignment"]
    num_sequences = len(final_alignments)

    # Clears existing widgets
    for widget in clustal_alignment.winfo_children():
        widget.destroy()

    # Style constants
    ROW_HEIGHT = 22
    ID_WIDTH = 25
    CHAR_WIDTH = 8
    MAX_VISIBLE_ROWS = 15
    ID_FONT = ("Courier New", 10, "bold")
    ALIGN_FONT = ("Courier New", 10)
    HEADER_FONT = ("Arial", 11, "bold")
    BG_COLOR = "#f5f5f5"
    HEADER_COLOR = "#3f51b5"
    TEXT_COLOR = "#333333"

    main_frame = ttk.Frame(clustal_alignment, padding=10)
    main_frame.pack(fill="both", expand=True)

    header_frame = ttk.Frame(main_frame)
    header_frame.pack(anchor="w", pady=(0, 10), fill="x")

    def make_label(text):
        return tk.Label(header_frame, text=text, font=HEADER_FONT, fg=HEADER_COLOR, bg=BG_COLOR)

    for label in [
        f"Num of Sequences: {num_sequences}",
        f"Scoring: {score}",
        f"Identity: {statistics.get('identity_percent')}%",
        f"Num of matches: {statistics.get('match')}",
        f"Num of mismatches: {statistics.get('mismatch')}",
        f"Num of gaps: {statistics.get('gap')}"
    ]:
        make_label(label).pack(side="left", padx=(0, 14))


    def save_result_in_clustal_format():
        """
           Displays and enables saving of the multiple sequence alignment in CLUSTAL format.

           Args:
               final_alignments (list): List of (name, alignment_string) tuples.
               score (int): Alignment score.
               statistics (dict): Contains identity_percent, match, mismatch, and gap counts.
               parameters_frame_dict (dict): Dictionary containing Entry widgets for alignment parameters.

           Functionality:
               - Clears previous CLUSTAL tab content.
               - Displays alignment metadata and statistics.
               - Allows saving results in CLUSTAL format to a text file.
               - Renders sequence alignments with scrollable UI components.
               - Synchronizes scrolling between ID panel and sequences.
           """
        file_path = filedialog.asksaveasfilename(
            defaultextension=".txt",
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")]
        )
        if file_path:
            max_name_length = max(len(name) for name, _ in final_alignments)
            with open(file_path, "w", encoding="utf-8") as file:

                file.write(f"Multiple Sequence Alignment in CLUSTAL Format\n")
                file.write(f"Statistics:\n")
                file.write(f"\tIdentity:   {statistics.get("identity_percent")}\n")
                file.write(f"\tScore:   {score}\n")
                file.write(f"\tNumber of Matches:   {statistics.get("match")}\n")
                file.write(f"\tNumber of MisMatches:   {statistics.get("mismatch")}\n")
                file.write(f"\tNumber of Gaps:   {statistics.get("gap")}\n\n")

                file.write(f"Matrix Scores:\n")
                file.write(f"\tMatch Score:   {parameters_frame_dict["match_score"].get()}\n")
                file.write(f"\tMismatch Score:   {parameters_frame_dict["mismatch_score"].get()}\n")
                file.write(f"\tGap Penalty:   {parameters_frame_dict["gap_score"].get()}\n\n")

                file.write(f"Scoring Result:\n")
                file.write(f"\tMatch:   {parameters_frame_dict["match"].get()}\n")
                file.write(f"\tSubstitution:   {parameters_frame_dict["substitution"].get()}\n")
                file.write(f"\tGap:   {parameters_frame_dict["gap"].get()}\n")

                file.write(f"Alignments:\n")
                for name, alignment in final_alignments:
                    formatted_line = f"{name.ljust(max_name_length)}  {''.join(alignment)}"
                    file.write(formatted_line + "\n")
    btn_frame = ttk.Frame(main_frame)
    btn_frame.pack(anchor="w", pady=(0, 20), fill="x")
    save_button = tk.Button(btn_frame, text="Save result in CLUSTAL Format", font=HEADER_FONT, fg=HEADER_COLOR,
            bg=BG_COLOR, command=save_result_in_clustal_format)
    save_button.pack(side="top", fill="both")

    # Creates container for horizontal scrollbar
    h_scroll_frame = ttk.Frame(main_frame)
    h_scroll_frame.pack(fill="x")

    # Horizontal scrollbar (top)
    h_scroll = ttk.Scrollbar(h_scroll_frame, orient="horizontal")
    h_scroll.pack(fill="x")

    # Create dual container for content
    content_frame = ttk.Frame(main_frame)
    content_frame.pack(fill="both", expand=True)

    # Vertical scrollbar (right side)
    v_scroll = ttk.Scrollbar(content_frame, orient="vertical")
    v_scroll.pack(side="right", fill="y")

    # Create canvas for scrolling
    canvas = tk.Canvas(
        content_frame,
        yscrollcommand=v_scroll.set,
        bg="white",
        highlightthickness=0,
        height=min(num_sequences, MAX_VISIBLE_ROWS) * ROW_HEIGHT,
    )
    canvas.pack(side="left", fill="both", expand=True)
    v_scroll.config(command=canvas.yview)

    # Create inner frame for IDs and sequences
    inner_frame = ttk.Frame(canvas)
    canvas.create_window((0, 0), window=inner_frame, anchor="nw")

    # Left panel for IDs
    id_frame = ttk.Frame(inner_frame, width=ID_WIDTH)
    id_frame.pack(side="left", fill="y")

    # Right panel for sequences
    seq_canvas = tk.Canvas(
        inner_frame,
        xscrollcommand=h_scroll.set,
        bg="white",
        highlightthickness=0
    )
    seq_canvas.pack(side="left", fill="both", expand=True)
    h_scroll.config(command=seq_canvas.xview)

    # Frame inside sequence canvas
    seq_frame = ttk.Frame(seq_canvas)
    seq_canvas.create_window((0, 0), window=seq_frame, anchor="nw")

    # Add IDs and sequences
    for seq_id, alignment in final_alignments:

        # Sequence ID
        id_label = tk.Label(
            id_frame,
            text=seq_id[:ID_WIDTH].ljust(ID_WIDTH),
            font=ID_FONT,
            fg=TEXT_COLOR,
            bg=BG_COLOR,
            anchor="w",
            padx=2,
        )
        id_label.pack(fill="x", pady=0)

        # Alignment sequence
        seq_label = tk.Label(
            seq_frame,
            text=alignment,
            font=ALIGN_FONT,
            fg=TEXT_COLOR,
            bg="white",
            anchor="w",
            padx=2
        )
        seq_label.pack(fill="x", pady=0, anchor="w")


    # Configure scrolling
    def configure_scrollregion(event=None):
        canvas.configure(scrollregion=canvas.bbox("all"))
        seq_canvas.configure(scrollregion=seq_canvas.bbox("all"))

        max_seq_len = max((len(al) for _, al in final_alignments), default=0)
        natural_width = max_seq_len * CHAR_WIDTH + 20  # Keep this
        available_width = clustal_alignment.winfo_width() - id_frame.winfo_width()  # Key change!
        seq_canvas.config(width=max(natural_width, available_width)) # Key change! Remove -50

    clustal_alignment.bind("<Configure>", configure_scrollregion)
    inner_frame.bind("<Configure>", configure_scrollregion)
    seq_frame.bind("<Configure>", configure_scrollregion)

    # Synchronized vertical scrolling
    def sync_yview(*args):
        canvas.yview(*args)
        seq_canvas.yview(*args)

    v_scroll.config(command=sync_yview)

    # Enhanced mouse/trackpad scrolling
    def _on_scroll(event):
        if seq_canvas.winfo_exists():
            if event.delta:
                seq_canvas.xview_scroll(int(-1 * (event.delta / 120)), "units")
            elif event.num == 6:
                seq_canvas.xview_scroll(1, "units")
            elif event.num == 7:
                seq_canvas.xview_scroll(-1, "units")
            else:
                sync_yview("scroll", int(-1 * (event.delta / 120)) if event.delta else -1, "units")

    seq_canvas.bind_all("<MouseWheel>", _on_scroll)
    seq_canvas.bind_all("<Button-4>", lambda e: _on_scroll(e))
    seq_canvas.bind_all("<Button-5>", lambda e: _on_scroll(e))
    seq_canvas.bind_all("<Shift-MouseWheel>", _on_scroll)

    # Final layout update
    inner_frame.update_idletasks()
    configure_scrollregion()

def print_result_in_alignment_viewer(final_alignments,score, statistics):
    """
       Display sequence alignments with color-coded nucleotides in a Tkinter viewer,
       including statistics and a button to save the alignment as a PNG image.

       Args:
           final_alignments (list of tuples): [(sequence_id, sequence_str), ...]
           score (int/float): Alignment score.
           statistics (dict): Alignment stats including keys:
               - 'identity_percent'
               - 'match'
               - 'mismatch'
               - 'gap'

       Behavior:
           - Clears and updates the Tkinter "AlignmentViewer" tab.
           - Shows alignment details and colored sequence.
           - Adds horizontal and vertical scrollbars.
           - Allows saving the alignment as a PNG file with colored bases.
       """
    viewer_alignment = tabs["AlignmentViewer"]
    num_sequences = len(final_alignments)

    with open("colors.json", "r") as f:
        color_data = json.load(f)
        AA_COLORS = {item["nucleotide"]: item["color"] for item in color_data}

    # Clear existing widgets
    for widget in viewer_alignment.winfo_children():
        widget.destroy()

    # Style constants
    ROW_HEIGHT = 22
    CHAR_WIDTH = 8
    ID_WIDTH = 25
    MAX_VISIBLE_ROWS = 15
    ID_FONT = ("Courier New", 10, "bold")
    ALIGN_FONT = ("Courier New", 12, "bold")
    HEADER_FONT = ("Arial", 11, "bold")
    BG_COLOR = "#f5f5f5"
    HEADER_COLOR = "#3f51b5"
    TEXT_COLOR = "#333333"

    # Main container
    main_frame = ttk.Frame(viewer_alignment, padding=10)
    main_frame.pack(fill="both", expand=True)

    # Header
    header_frame = ttk.Frame(main_frame)
    header_frame.pack(anchor="w", pady=(0, 10), fill="x")

    def make_label(text):
        return tk.Label(header_frame, text=text, font=HEADER_FONT, fg=HEADER_COLOR, bg=BG_COLOR)

    for label in [
        f"Num of Sequences: {num_sequences}",
        f"Scoring: {score}",
        f"Identity: {statistics.get('identity_percent')}%",
        f"Num of matches: {statistics.get('match')}",
        f"Num of mismatches: {statistics.get('mismatch')}",
        f"Num of gaps: {statistics.get('gap')}"
    ]:
        make_label(label).pack(side="left", padx=(0, 14))


    def save_result_as_png():
        with open("colors.json", "r") as f:
            raw_colors = json.load(f)
            color_map = {entry["nucleotide"].upper(): entry["color"] for entry in raw_colors}

        cell_size = 30
        name_width = 150
        padding = 5
        font = ImageFont.load_default()
        font_name = ImageFont.truetype("arial.ttf", 15)
        font_score = ImageFont.truetype("arial.ttf", 13)

        num_rows = len(final_alignments)
        num_cols = max(len(seq) for _, seq in final_alignments)

        width = name_width + num_cols * cell_size
        height = (num_rows + 1) * cell_size  # extra row for the score

        img = Image.new("RGB", (width, height), "white")
        draw = ImageDraw.Draw(img)

        num_sequences = len(final_alignments)

        padding_scoring_x = padding
        padding_scoring_y = padding

        for text in [
            f"Num of Sequences: {num_sequences}",
            f"Scoring: {score}",
            f"Identity: {statistics.get('identity_percent')}%",
            f"Matches: {statistics.get('match')}",
            f"Mismatches: {statistics.get('mismatch')}",
            f"Gaps: {statistics.get('gap')}",
        ]:
            draw.text((padding_scoring_x, padding_scoring_y), text, fill="black", font=font_score)
            text_width = draw.textlength(text, font=font_score)
            padding_scoring_x += text_width + 20  # adds space between blocks

        for row_idx, (name, sequence) in enumerate(final_alignments):
            y = (row_idx + 1) * cell_size  # shifts down for score

            # Draws the sequence name
            draw.text((5, y + padding), name[:15], fill="black", font=font_name)

            for col_idx, char in enumerate(sequence):
                x = name_width + col_idx * cell_size
                color = color_map.get(char.upper(), "#FFFFFF")  # fallback to white

                # Draws colored cell
                draw.rectangle([x, y, x + cell_size, y + cell_size], fill=color, outline="black")

                # Draws centered character
                bbox = draw.textbbox((0, 0), char, font=font)
                w = bbox[2] - bbox[0]
                h = bbox[3] - bbox[1]
                draw.text(
                    (x + (cell_size - w) / 2, y + (cell_size - h) / 2),
                    char, fill="black", font=font
                )
        img.save("alignment.png")
        # Opens a save dialog to choose destination and filename
        file_path = filedialog.asksaveasfilename(
            defaultextension=".png",
            filetypes=[("PNG files", "*.png")])

        if file_path:
            img.save(file_path)

        # Header
    btn_frame = ttk.Frame(main_frame)
    btn_frame.pack(anchor="w", pady=(0, 20), fill="x")
    save_button = tk.Button(btn_frame, text="Save result as .png", font=HEADER_FONT, fg=HEADER_COLOR,
            bg=BG_COLOR, command=save_result_as_png)
    save_button.pack(side="top", fill="both")

    # Creates container for horizontal scrollbar
    h_scroll_frame = ttk.Frame(main_frame)
    h_scroll_frame.pack(fill="x")

    # Horizontal scrollbar (top)
    h_scroll = ttk.Scrollbar(h_scroll_frame, orient="horizontal")
    h_scroll.pack(fill="x")

    # Creates dual container for content
    content_frame = ttk.Frame(main_frame)
    content_frame.pack(fill="both", expand=True)

    # Vertical scrollbar (right side)
    v_scroll = ttk.Scrollbar(content_frame, orient="vertical")
    v_scroll.pack(side="right", fill="y")

    # Creates canvas for scrolling
    canvas = tk.Canvas(
        content_frame,
        yscrollcommand=v_scroll.set,
        bg="white",
        highlightthickness=0,
        height=min(num_sequences, MAX_VISIBLE_ROWS) * ROW_HEIGHT,
    )
    canvas.pack(side="left", fill="both", expand=True)
    v_scroll.config(command=canvas.yview)

    # Create inner frame for IDs and sequences
    inner_frame = ttk.Frame(canvas)
    canvas.create_window((0, 0), window=inner_frame, anchor="nw")

    # Left panel for IDs
    id_frame = ttk.Frame(inner_frame, width=ID_WIDTH)
    id_frame.pack(side="left", fill="y")

    # Right panel for sequences
    seq_canvas = tk.Canvas(
        inner_frame,
        xscrollcommand=h_scroll.set,
        bg="white",
        highlightthickness=0
    )
    seq_canvas.pack(side="left", fill="both", expand=True)
    h_scroll.config(command=seq_canvas.xview)

    # Frame inside sequence canvas
    seq_frame = ttk.Frame(seq_canvas)
    seq_canvas.create_window((0, 0), window=seq_frame, anchor="nw")

    # Adds IDs and sequences
    for seq_id, alignment in final_alignments:
        # Sequence ID
        id_label = tk.Label(
            id_frame,
            text=seq_id[:ID_WIDTH].ljust(ID_WIDTH),
            font=ID_FONT,
            fg=TEXT_COLOR,
            bg=BG_COLOR,
            anchor="w",
            padx=2,
        )
        id_label.pack(fill="x", pady=0)

        # Colored character-by-character alignment display
        char_row = ttk.Frame(seq_frame)
        char_row.pack(anchor="w", pady=0)

        for char in alignment:
            bg = AA_COLORS.get(char.upper(), "white")
            char_label = tk.Label(
                char_row,
                text=char,
                font=ALIGN_FONT,
                fg="black",
                bg=bg,
                width=2,
                padx=1,
                pady=1,
                borderwidth=1,
                relief="solid"
            )
            char_label.pack(side="left")

    # Configure scrolling
    def configure_scrollregion(event=None):
        canvas.configure(scrollregion=canvas.bbox("all"))
        seq_canvas.configure(scrollregion=seq_canvas.bbox("all"))

        max_seq_len = max((len(al) for _, al in final_alignments), default=0)
        natural_width = max_seq_len * CHAR_WIDTH + 20  # Keep this
        available_width = viewer_alignment.winfo_width() - id_frame.winfo_width()  # Key change!
        seq_canvas.config(width=max(natural_width, available_width)) # Key change! Remove -50



    viewer_alignment.bind("<Configure>", configure_scrollregion)
    inner_frame.bind("<Configure>", configure_scrollregion)
    seq_frame.bind("<Configure>", configure_scrollregion)

    # Synchronized vertical scrolling
    def sync_yview(*args):
        canvas.yview(*args)
        seq_canvas.yview(*args)

    v_scroll.config(command=sync_yview)

    # Enhanced mouse/trackpad scrolling
    def _on_scroll(event):
        if seq_canvas.winfo_exists():
            if event.delta:
                seq_canvas.xview_scroll(int(-1 * (event.delta / 120)), "units")
            elif event.num == 6:
                seq_canvas.xview_scroll(1, "units")
            elif event.num == 7:
                seq_canvas.xview_scroll(-1, "units")
            else:
                sync_yview("scroll", int(-1 * (event.delta / 120)) if event.delta else -1, "units")

    seq_canvas.bind_all("<MouseWheel>", _on_scroll)
    seq_canvas.bind_all("<Button-4>", lambda e: _on_scroll(e))
    seq_canvas.bind_all("<Button-5>", lambda e: _on_scroll(e))
    seq_canvas.bind_all("<Shift-MouseWheel>", _on_scroll)

    # Final layout update
    inner_frame.update_idletasks()
    configure_scrollregion()

def print_result_in_fasta_format(final_alignments,score,statistics, parameters_frame_dict):
    """
        Display sequence alignments in FASTA format inside a Tkinter tab,
        showing alignment stats and providing a save-to-file option.

        Args:
            final_alignments (list of tuples): [(sequence_id, sequence_str), ...]
            score (int/float): Alignment score.
            statistics (dict): Alignment statistics with keys:
                - 'identity_percent', 'match', 'mismatch', 'gap'
            parameters_frame_dict (dict): GUI elements holding scoring parameters.

        Behavior:
            - Clears and updates the "FASTA Alignment" tab.
            - Shows sequence headers and sequences formatted with line breaks.
            - Displays alignment stats and scoring parameters.
            - Allows saving the alignment and stats as a text file in FASTA format.
        """
    fasta_alignment = tabs["FASTA Alignment"]
    num_sequences = len(final_alignments)

    # Clear existing widgets
    for widget in fasta_alignment.winfo_children():
        widget.destroy()

    # Style constants
    ROW_HEIGHT = 22
    CHAR_WIDTH = 8
    ID_WIDTH = 25
    MAX_VISIBLE_ROWS = 15
    ID_FONT = ("Courier New", 10, "bold")
    SEQUENCE_FONT = ("Courier New", 10)
    HEADER_FONT = ("Arial", 11, "bold")
    BG_COLOR = "#f5f5f5"
    HEADER_COLOR = "#3f51b5"
    TEXT_COLOR = "#333333"
    RIGHT_PADDING = 18  # Extra space to prevent hiding behind scrollbar
    LINE_LENGTH = 55  # Standard FASTA line length

    # Main container
    main_frame = ttk.Frame(fasta_alignment, padding=10)
    main_frame.pack(fill="both", expand=True)

    # Header
    header_frame = ttk.Frame(main_frame)
    header_frame.pack(anchor="w", pady=(0, 10), fill="x")

    def make_label(text):
        return tk.Label(header_frame, text=text, font=HEADER_FONT, fg=HEADER_COLOR, bg=BG_COLOR)

    for label in [
        f"Num of Sequences: {num_sequences}",
        f"Scoring: {score}",
        f"Identity: {statistics.get('identity_percent')}%",
        f"Num of matches: {statistics.get('match')}",
        f"Num of mismatches: {statistics.get('mismatch')}",
        f"Num of gaps: {statistics.get('gap')}"
    ]:
        make_label(label).pack(side="left", padx=(0, 14))

    def save_result_in_fasta_format():
        # Saves data to the selected file path
        file_path = filedialog.asksaveasfilename(
            defaultextension=".txt",
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")]
        )
        if file_path:
            # First, find the max length of all names
            max_name_length = max(len(name) for name, _ in final_alignments)
            # Now write the file
            with open(file_path, "w", encoding="utf-8") as file:

                file.write(f"Multiple Sequence Alignment in FASTA Format\n")
                file.write(f"Statistics:\n")
                file.write(f"\tIdentity:   {statistics.get("identity_percent")}\n")
                file.write(f"\tScore:   {score}\n")
                file.write(f"\tNumber of Matches:   {statistics.get("match")}\n")
                file.write(f"\tNumber of MisMatches:   {statistics.get("mismatch")}\n")
                file.write(f"\tNumber of Gaps:   {statistics.get("gap")}\n")


                file.write(f"Matrix Scores:\n")
                file.write(f"\tMatch Score:   {parameters_frame_dict["match_score"].get()}\n")
                file.write(f"\tMismatch Score:   {parameters_frame_dict["mismatch_score"].get()}\n")
                file.write(f"\tGap Penalty:   {parameters_frame_dict["gap_score"].get()}\n")

                file.write(f"Scoring Result:\n")
                file.write(f"\tMatch:   {parameters_frame_dict["match"].get()}\n")
                file.write(f"\tSubstitution:   {parameters_frame_dict["substitution"].get()}\n")
                file.write(f"\tGap:   {parameters_frame_dict["gap"].get()}\n")


                file.write(f"Alignments:\n")
                for name, alignment in final_alignments:
                    file.write(f">{name}\n")
                    formatted = (f"{''.join(alignment)}\n")
                    file.write(formatted)
                #file.write(df.to_string(index=False, header=False))
        # Header
    btn_frame = ttk.Frame(main_frame)
    btn_frame.pack(anchor="w", pady=(0, 20), fill="x")
    save_button = tk.Button(btn_frame, text="Save result in FASTA Format", font=HEADER_FONT, fg=HEADER_COLOR,
            bg=BG_COLOR, command=save_result_in_fasta_format)
    save_button.pack(side="top", fill="both")

    # Create container for content
    content_frame = ttk.Frame(main_frame)
    content_frame.pack(fill="both", expand=True)

    # Vertical scrollbar
    v_scroll = ttk.Scrollbar(content_frame, orient="vertical")
    v_scroll.pack(side="right", fill="y")

    # Create canvas for scrolling
    canvas = tk.Canvas(
        content_frame,
        yscrollcommand=v_scroll.set,
        bg="white",
        highlightthickness=0,
        height=min(num_sequences * 3, MAX_VISIBLE_ROWS * 3) * ROW_HEIGHT
    )
    canvas.pack(side="left", fill="both", expand=True, padx=(0, RIGHT_PADDING))
    v_scroll.config(command=canvas.yview)

    # Create inner frame for FASTA content
    inner_frame = ttk.Frame(canvas)
    canvas.create_window((0, 0), window=inner_frame, anchor="nw")

    # Calculate maximum width needed
    def calculate_max_width():
        max_header = max(len(f">{seq_id}") for seq_id, _ in final_alignments) if final_alignments else 0
        max_sequence = LINE_LENGTH  # We're controlling line length to 60 chars
        return max(max_header, max_sequence) * CHAR_WIDTH + RIGHT_PADDING

    max_width = calculate_max_width()

    # Configure inner frame width
    inner_frame.config(width=max_width)

    # Add FASTA formatted sequences
    for seq_id, alignment in final_alignments:
        # FASTA header line
        header_label = tk.Label(
            inner_frame,
            text=f">{seq_id}",
            font=ID_FONT,
            fg=HEADER_COLOR,
            bg=BG_COLOR,
            anchor="w",
            justify="left",
            width= 110  # Match sequence line length for alignment
        )
        header_label.pack(fill="x", anchor="w")

        # Sequence lines (split into standard FASTA chunks)
        for i in range(0, len(alignment), LINE_LENGTH):
            chunk = alignment[i:i+LINE_LENGTH]
            seq_label = tk.Label(
                inner_frame,
                text=chunk,
                font=SEQUENCE_FONT,
                fg=TEXT_COLOR,
                bg="white",
                anchor="w",
                justify="left",
                width=LINE_LENGTH  # Fixed width for all lines
            )
            seq_label.pack(fill="x", anchor="w")

        # Add empty line between sequences
        spacer = tk.Label(inner_frame, text="", bg="white", height=1)
        spacer.pack(fill="x")

    # Configure scrolling with proper padding
    def configure_scrollregion(event):
        # Get the required size including all content
        bbox = canvas.bbox("all")
        if bbox:
            # Add right padding to ensure no content is hidden
            canvas.configure(scrollregion=(0, 0, bbox[2] + RIGHT_PADDING, bbox[3]))
            # Set canvas width to show full content
            canvas.config(width=min(bbox[2] + RIGHT_PADDING, main_frame.winfo_width()))

    inner_frame.bind("<Configure>", configure_scrollregion)
    main_frame.bind("<Configure>", lambda e: configure_scrollregion(None))

    # Mouse wheel scrolling
    def _on_mousewheel(event):
        canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")

    canvas.bind_all("<MouseWheel>", _on_mousewheel)

    # Initial configuration
    inner_frame.update_idletasks()
    canvas.configure(scrollregion=canvas.bbox("all"))


def load_fasta_or_folder(input_frame_dict):
    """
        Opens a file dialog to select a single .fasta file or a folder containing .fasta files,
        loads their contents into a global list `fasta_contents`, and updates UI message.

        Args:
            input_frame_dict (dict): Contains UI elements, including 'load_file_message' Label for status updates.

        Behavior:
            - If user cancels file selection, prompts to select a folder.
            - Loads .fasta file(s) contents into `fasta_contents`.
            - Shows error if invalid file or no FASTA files in folder.
        """
    selected_path = filedialog.askopenfilename(title="Select a .fasta file or cancel to select folder")

    if not selected_path:  # user canceled file selection, ask folder instead
        selected_path = filedialog.askdirectory(title="Select .fasta file or folder containing .fasta files")
        if not selected_path:
            return  # user canceled folder selection too, do nothing

    global fasta_contents
    if os.path.isfile(selected_path):
        if selected_path.lower().endswith(".fasta"):
            with open(selected_path, 'r') as f:
                fasta_contents.append(f.read())
            input_frame_dict['load_file_message'].config(text=f"Loaded FASTA file: {selected_path}")
        else:
            messagebox.showerror("Invalid file", "Selected file is not a .fasta file.")
            return
    elif os.path.isdir(selected_path):
        # It's a folder - load all .fasta files inside
        fasta_files = [f for f in os.listdir(selected_path) if f.lower().endswith(".fasta")]
        if not fasta_files:
            messagebox.showwarning("No FASTA files", "No .fasta files found in the selected folder.")
            return
        for fasta_file in fasta_files:
            full_path = os.path.join(selected_path, fasta_file)
            with open(full_path, 'r') as f:
                fasta_contents.append(f.read())


        """input_frame_dict['load_file_message'].config(text=f"Loaded folder with FASTA files: {selected_path}")
        input_frame_dict['load_file_message'].pack(fill="x")"""
    else:
        messagebox.showerror("Invalid selection", "Selected path is not a file or folder.")
        return


def extract_names_and_sequences(fasta_content):
    """
    Extracts sequence names (headers) and their sequences from FASTA content.

    :param fasta_content: string, the whole content of a FASTA file
    :return: list of tuples (header_name, sequence)
    """
    sequences = {}
    current_name = None
    current_seq_lines = []

    for line in fasta_content.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            # Save previous sequence if any
            if current_name:
                sequences[current_name] = "".join(current_seq_lines)
            # Extract name (header line without '>')
            current_name = line[1:].strip()
            current_seq_lines = []
        else:
            # Validate sequence line - only amino acid letters
            if re.fullmatch(r'[ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy]+', line):
                current_seq_lines.append(line)
            else:
                messagebox.showerror("Invalid sequence", f"The sequence {current_name} contains invalid characters. Please correct it.")
                raise Exception("Invalid sequence")

    # Save the last sequence
    if current_name:
        sequences[current_name] = "".join(current_seq_lines)

    # Convert to list of tuples and return
    return list(sequences.items())


def main():

    # Input tab
    input_frame = tabs["Input"]

    text_container = tk.Frame(input_frame, width=1000, height=380)
    text_container.pack(pady=(10, 5))
    text_container.pack_propagate(False)  # Prevent auto-resizing

    scrollbar_y = tk.Scrollbar(text_container, orient="vertical")
    scrollbar_y.pack(side="right", fill="y")

    scrollbar_x = tk.Scrollbar(text_container, orient="horizontal")
    scrollbar_x.pack(side="bottom", fill="x")

    text_field = tk.Text(text_container, wrap="none",
                         height=15, font=("Courier New", 11), fg="grey",
                         yscrollcommand=scrollbar_y.set,
                         xscrollcommand=scrollbar_x.set)
    text_field.pack(side="left", fill="both", expand=True)


    scrollbar_y.config(command=text_field.yview)
    scrollbar_x.config(command=text_field.xview)

    bind_mousewheel(text_field, text_field)

    setup_placeholder(text_field)

    btn_frame_input = ttk.Frame(input_frame)
    btn_frame_input.pack(side="bottom", anchor="e", pady=10, padx=10)

    button_style = {
        "font": ("Segoe UI", 10, "bold"),
        "bd": 0,
        "relief": "flat",
        "padx": 20,
        "pady": 10,
        "cursor": "hand2"
    }

    submit_btn = tk.Button(
        btn_frame_input,
        text="Submit",
        bg="#437871",
        fg="white",
        activebackground="#256f47",
        activeforeground="white",
        **button_style
    )
    submit_btn.pack(side="right", padx=(10, 10))

    reset_btn = tk.Button(
        btn_frame_input,
        text="Reset",
        bg="#cc4c4c",
        fg="white",
        activebackground="#a23a3a",
        activeforeground="white",
        **button_style
    )
    reset_btn.pack(side="right", padx=(0, 10))

    load_file_message = tk.Label(input_frame, text="", font=("Arial", 12))
    load_file_message.pack(side="right", padx=(0, 10))

    input_frame_dict = {
        'text_field': text_field,
        'submit_btn': submit_btn,
        'reset_btn': reset_btn,
        'load_file_message': load_file_message,
    }

    submit_btn.config(command=lambda: on_submit_btn_click(input_frame_dict, parameters_frame_dict))
    reset_btn.config(command=lambda: on_reset_btn_click(input_frame_dict, parameters_frame_dict))

    def paste_example():
        with open("example_sequences", "r") as f:
            example_sequences = f.read()
        text_field.delete("1.0", "end")
        text_field.insert("1.0", example_sequences)
        text_field.config(fg="black")

    btn_frame = tk.Frame(input_frame, bg="#f4f4f4")
    btn_frame.pack(anchor="w", pady=(0, 5))

    paste_btn = tk.Button(btn_frame, command=paste_example, text="Paste Example",
                          fg="#3e8074", bg="#f4f4f4", bd=0,
                          font=("Arial", 12), cursor="hand2")
    paste_btn.pack(side="left", padx=(0, 20))

    upload_btn = tk.Button(btn_frame, text="Upload File",
                           fg="#3e8074", bg="#f4f4f4", bd=0,
                           font=("Arial", 12, ), cursor="hand2")
    upload_btn.pack(side="left")

    upload_btn.config(command=lambda: load_fasta_or_folder(input_frame_dict))

    # Parameters tab
    param_frame = tabs["Parameters"]

    # Parameters Frame
    parameters_frame = ttk.Frame(param_frame)
    parameters_frame.pack(anchor="w", pady=20)

    # --- Matrix Score row ---
    matrix_label = ttk.Label(parameters_frame, text="Matrix Score", font=("Arial", 12, "bold"))
    matrix_label.grid(row=0, column=0, sticky="w", pady=(0, 2))

    # Labels above matrix entries
    ttk.Label(parameters_frame, text="Match").grid(row=1, column=0)
    ttk.Label(parameters_frame, text="Mismatch").grid(row=1, column=1)
    ttk.Label(parameters_frame, text="Gap").grid(row=1, column=2)

    # Matrix Score entries
    match_score_entry = tk.Entry(parameters_frame, width=10, font=("Arial", 12))
    match_score_entry.grid(row=2, column=0, padx=5)
    match_score_entry.insert(0, "1")

    mismatch_score_entry = tk.Entry(parameters_frame, width=10, font=("Arial", 12))
    mismatch_score_entry.grid(row=2, column=1, padx=5)
    mismatch_score_entry.insert(0, "-1")

    gap_score_entry = tk.Entry(parameters_frame, width=10, font=("Arial", 12))
    gap_score_entry.grid(row=2, column=2, padx=5)
    gap_score_entry.insert(0, "-2")

    # Scoring Result row
    result_label = ttk.Label(parameters_frame, text="Scoring Result", font=("Arial", 12, "bold"))
    result_label.grid(row=3, column=0, sticky="w", pady=(10, 2))

    # Labels above scoring result entries
    ttk.Label(parameters_frame, text="Match").grid(row=4, column=0)
    ttk.Label(parameters_frame, text="Substitution").grid(row=4, column=1)
    ttk.Label(parameters_frame, text="Gap").grid(row=4, column=2)

    # Scoring Result entries
    match_entry = tk.Entry(parameters_frame, width=10, font=("Arial", 12))
    match_entry.grid(row=5, column=0, padx=5)
    match_entry.insert(0, "1")

    substitution_entry = tk.Entry(parameters_frame, width=10, font=("Arial", 12))
    substitution_entry.grid(row=5, column=1, padx=5)
    substitution_entry.insert(0, "-1")

    gap_entry = tk.Entry(parameters_frame, width=10, font=("Arial", 12))
    gap_entry.grid(row=5, column=2, padx=5)
    gap_entry.insert(0, "-2")

    # Frame for buttons inside the same 'param_frame'
    btn_frame_param = ttk.Frame(param_frame)
    btn_frame_param.pack(side="bottom", anchor="e", pady=10, padx=10)

    button_style = {
        "font": ("Segoe UI", 10, "bold"),
        "bd": 0,
        "relief": "flat",
        "padx": 20,
        "pady": 10,
        "cursor": "hand2"
    }

    resubmit_btn = tk.Button(
        btn_frame_param,
        text="Submit",
        bg="#437871",
        fg="white",
        activebackground="#256f47",
        activeforeground="white",
        **button_style
    )
    resubmit_btn.pack(side="right", padx=(10, 10))

    reset_btn2 = tk.Button(
        btn_frame_param,
        text="Reset",
        bg="#cc4c4c",
        fg="white",
        activebackground="#a23a3a",
        activeforeground="white",
        **button_style
    )
    reset_btn2.pack(side="right", padx=(0, 10))

    parameters_frame_dict = {
        'reset_btn2': reset_btn2,
        'resubmit_btn': resubmit_btn,
        'match_score': match_score_entry,
        'mismatch_score': mismatch_score_entry,
        'gap_score': gap_score_entry,
        'match' : match_entry,
        'substitution' : substitution_entry,
        'gap' : gap_entry
    }

    resubmit_btn.config(command=lambda: on_submit_btn_click(input_frame_dict, parameters_frame_dict))
    reset_btn2.config(command=lambda: on_reset_btn_click(input_frame_dict, parameters_frame_dict))

    root.mainloop()
import os
import sys
import re

def extract_strings_from_binary(file_path, min_length=4):
    with open(file_path, 'rb') as f:
        data = f.read()
    # Find all printable sequences
    result = ""
    for match in re.finditer(b'[a-zA-Z0-9\s\.\,\:\;\-\_\(\)\{\}\[\]\;\?\!]{' + str(min_length).encode() + b',}', data):
        try:
            s = match.group().decode('utf-8', errors='ignore')
            result += s + "\n"
        except:
            pass
    return result

def extract_pdf(file_path):
    try:
        from pypdf import PdfReader
    except ImportError:
        try:
            import PyPDF2 as PdfReader
        except ImportError:
            return "ERROR: No PDF library found (pypdf or PyPDF2)."
    
    try:
        reader = PdfReader(file_path)
        text = ""
        for page in reader.pages:
            text += page.extract_text() + "\n"
        return text
    except Exception as e:
        return f"ERROR reading PDF: {e}"

def main():
    base_dir = r"c:\Users\Goksun\MATLAB\Projects\fhss_senior_project"
    pdf_path = os.path.join(base_dir, "EE491_492_Midterm_Report_GoksunGuneyKucuk.pdf")
    doc_path = os.path.join(base_dir, "EE491-EE492-final_report (1).doc")

    print("--- EXTRACTING PDF ---")
    pdf_text = extract_pdf(pdf_path)
    print(pdf_text[:5000]) # First 5000 chars roughly
    print("\n... [truncated] ...\n")
    
    print("--- EXTRACTING DOC STRINGS ---")
    doc_text = extract_strings_from_binary(doc_path)
    # The doc likely has a lot of noise, so print a chunk to see if we get the template headers
    print(doc_text[:5000])

if __name__ == "__main__":
    main()

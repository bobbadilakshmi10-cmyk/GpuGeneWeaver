
# Import Flask and required modules
from flask import Flask, render_template, request, jsonify

# Import BioPython to read FASTA DNA sequence files
from Bio import SeqIO

# Import StringIO to read uploaded text as a file
from io import StringIO


# Create Flask application
app = Flask(__name__)


# Home page route
@app.route("/")
def home():
    # Display the index.html webpage
    return render_template("index.html")


# API health check route
@app.route("/api/health")
def health():

    # Return application status in JSON format
    return jsonify({
        "status": "success",
        "message": "GeneWeaver is running"
    })


# DNA analysis route
# POST method is used to upload FASTA files and target DNA
@app.route("/api/analyze", methods=["POST"])
def analyze():

    # Get the uploaded FASTA file from the request
    fasta_file = request.files.get("fasta_file")

    # Get target DNA sequence and convert it to uppercase
    target = request.form.get("target", "").upper()

    # Get maximum allowed mismatches
    # Default mismatch value is 3
    max_mismatches = int(request.form.get("max_mismatches", 3))


    # Check whether FASTA file and target DNA are provided
    if not fasta_file or not target:
        return jsonify({
            "error": "Upload FASTA file and enter target DNA"
        }), 400


    try:
        # Read uploaded FASTA file content
        # utf-8-sig also handles a possible UTF-8 BOM
        content = fasta_file.read().decode("utf-8-sig")

        # Parse FASTA content using BioPython
        # StringIO converts text into a file-like object
        records = list(SeqIO.parse(StringIO(content), "fasta"))

        # Create an empty list to store matching results
        results = []


        # Process each DNA record in the FASTA file
        for record in records:

            # Convert DNA sequence into uppercase
            sequence = str(record.seq).upper()


            # Slide through the DNA sequence
            # Compare chunks having the same length as target
            for i in range(len(sequence) - len(target) + 1):

                # Extract a DNA chunk from the current position
                chunk = sequence[i:i + len(target)]


                # Count mismatches between target and DNA chunk
                # zip() compares characters position by position
                mismatches = sum(
                    1 for a, b in zip(chunk, target) if a != b
                )


                # Check whether mismatches are within the allowed limit
                if mismatches <= max_mismatches:

                    # Store matching information in results list
                    results.append({
                        "sequence_id": record.id,  # FASTA sequence ID
                        "position": i + 1,          # 1-based DNA position
                        "match": chunk,             # Matching DNA chunk
                        "mismatches": mismatches    # Number of mismatches
                    })


        # Return successful analysis results in JSON format
        return jsonify({
            "status": "success",
            "total_matches": len(results),
            "results": results
        })


    # Handle errors during file reading or DNA analysis
    except Exception as e:

        # Return error message with HTTP 400 status
        return jsonify({"error": str(e)}), 400


# Main entry point of the Python application
if __name__ == "__main__":

    # Start Flask development server
    # debug=True enables automatic reload and debugging
    app.run(debug=True)

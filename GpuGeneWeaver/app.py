
from flask import Flask, render_template, request, jsonify
from Bio import SeqIO
from io import StringIO

app = Flask(__name__)

@app.route("/")
def home():
    return render_template("index.html")

@app.route("/api/health")
def health():
    return jsonify({
        "status": "success",
        "message": "GeneWeaver is running"
    })

@app.route("/api/analyze", methods=["POST"])
def analyze():
    fasta_file = request.files.get("fasta_file")
    target = request.form.get("target", "").upper()
    max_mismatches = int(request.form.get("max_mismatches", 3))

    if not fasta_file or not target:
        return jsonify({"error": "Upload FASTA file and enter target DNA"}), 400

    try:
        content = fasta_file.read().decode("utf-8-sig")
        records = list(SeqIO.parse(StringIO(content), "fasta"))

        results = []

        for record in records:
            sequence = str(record.seq).upper()

            for i in range(len(sequence) - len(target) + 1):
                chunk = sequence[i:i + len(target)]

                mismatches = sum(
                    1 for a, b in zip(chunk, target) if a != b
                )

                if mismatches <= max_mismatches:
                    results.append({
                        "sequence_id": record.id,
                        "position": i + 1,
                        "match": chunk,
                        "mismatches": mismatches
                    })

        return jsonify({
            "status": "success",
            "total_matches": len(results),
            "results": results
        })

    except Exception as e:
        return jsonify({"error": str(e)}), 400


if __name__ == "__main__":
    app.run(debug=True)

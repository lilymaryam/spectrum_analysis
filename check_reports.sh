for file in ./reports/*; do
  if [ -f "$file" ]; then
    line_count=$(wc -l < "$file")
    second_col=$(awk 'NR==2{print $1}' "$file")

    echo "=== $file ==="

    if [ "$line_count" -gt 2 ]; then
      echo "  ✓ Has more than 2 lines ($line_count)"
      pass
    #else
    #  echo "  ✗ Only $line_count line(s)"
    fi

    #if [ "$second_col" = "node_1" ]; then
    #  echo "  ✓ Line 2, col 1 is 'node_1'"
    if [ "$second_col" != "node_1" ]; then
      echo "  ✗ Line 2, col 1 is '$second_col' (expected 'node_1')"
    fi
  fi
done
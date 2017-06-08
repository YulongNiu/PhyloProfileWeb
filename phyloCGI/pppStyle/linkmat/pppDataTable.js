d3.text("./predicted_linakges.csv", function(text) {

    var linksArray = d3.csvParseRows(text);
    linksArray.shift();

    var table = $(document).ready(function() {
        $('#linksmat').DataTable( {
            data: linksArray,
            columns: [
                { title: "From" },
                { title: "FromAnno" },
                { title: "To" },
                { title: "ToAnno" },
                { title: "Jaccard" },
                { title: "Cor" }
            ]
        } );
    } );

    $('#linksmat tbody')
        .on( 'mouseenter', 'td', function () {
            var colIdx = table.cell(this).index().column;

            $( table.cells().nodes() ).removeClass( 'highlight' );
            $( table.column( colIdx ).nodes() ).addClass( 'highlight' );
        } );
});

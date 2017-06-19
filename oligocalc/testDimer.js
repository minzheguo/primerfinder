function stringToArray(theString)
{	
	theArray=new Array(theString.length);
	for( var i=0; i<theString.length; i++) 
		theArray[i]=theString.charAt(i); 
	return theArray;
}


function MakeComplement(theSequence, isDNA) {
	var returnString="";
	var i;
	var temp;
	for( i=theSequence.length-1; i>=0; i--) {
		temp=theSequence.charAt(i);
		switch (temp) {
			case "A" :
				if (isDNA) {
					temp="T";
				} else {
					temp="U";
				}
				break;
			case "T" :
				temp="A";
				break;
			case "U" :
				temp="A";
				break;
			case "G" :
				temp="C";
				break;
			case "C" :
				temp="G";
				break;
			case "M" :
				temp="K";
				break;
			case "K" :
				temp="M";
				break;
			case "R" :
				temp="Y";
				break;
			case "Y" :
				temp="R";
				break;
			case "W" :
				temp="W";
				break;
			case "S" :
				temp="S";
				break;
			case "V" :
				temp="B";
				break;
			case "B" :
				temp="V";
				break;
			case "H" :
				temp="D";
				break;
			case "D" :
				temp="H";
				break;
			default : break;
		}
		returnString=returnString+temp;
	}
	return returnString;
}


function makeMatrix(matLength) 
{
	var theMatrix = new Array(matLength);
	for(var i=0; i<matLength; i++) {
		// increment column
		theMatrix[i]=new Array(matLength);
	}
	return theMatrix;
}

function fillMatchMatrix(cols, rows, mat, broadMatch=false)
{
	var d=cols.length;
	if (d<4) return;
	if (broadMatch) {
		alert("broadMatch")
		// Do the degenerate thing!
		for(i=0; i<d; i++) {
			// increment column
			for(var j=0; j<d; j++) {
				// increment row
				if(isBaseEqual( cols[i], rows[j]) ) {
					mat[i][j]=1;
					if(i>0 && j>0)
					{
						mat[i][j]+=mat[i-1][j-1]; //(increment diagonal values)
					}
				} else {	
					mat[i][j]=0;
					if(i>1 && j>1 ) {
						if(mat[i-1][j-1]>mat[i-2][j-2] &&mat[i-1][j-1]>1 && i<d-1 && j<d-1) {
						// allow one base mismatch only if there are at least 2 matched base on 5' and at least 1 matched base on 3'
							mat[i][j]=mat[i-1][j-1]; 
						} else if(i<d-1 && j<d-1) {
							mat[i-1][j-1]=0;
						}
					}
				}
			}
		}
	} else {
		for (i=0; i<=1; i++) {
			// increment column
			for (var j=0; j<2; j++) {
				// increment row
				if (cols[i] == rows[j]) {
					mat[i][j]=1;
					if (i && j)
						mat[i][j]+=mat[i-1][j-1]; //(increment diagonal values)
				} else {	
					mat[i][j]=0;
				}
			}
		}
		for (i=2; i<d-1; i++) {
			// increment column
			for (var j=2; j<d-1; j++) {
				// increment row
				if (cols[i] == rows[j]) {
					mat[i][j]=mat[i-1][j-1]+1; //(increment diagonal values)
				} else {	
					mat[i][j]=0;
					if (mat[i-1][j-1]>1 && cols[i+1] == rows[j+1]) {
					// allow one base mismatch only if there are at least 2 matched base on 5' and at least 1 matched base on 3'
						mat[i][j]=mat[i-1][j-1]; 
					}
				}
			}
		}
		i=d-1;
		j=i;
		// increment column
		// increment row
		if (cols[i] == rows[j]) {
			mat[i][j]=1;
			mat[i][j]+=mat[i-1][j-1]; //(increment diagonal values)
		} else {	
			mat[i][j]=0;
		}
	}
}

function makeAlignedArray(mat, minLen, maxMisMatch)
{
// assumes an orthogonal matrix
/* theAlignedArray is a bit strange in the second dimension. Assume it is a length 5 array called 'theResults'
	theResults[0] == start index
	theResults[1] == start matching index in reverse complement seq
	theResults[2] == end index of aligned bases (inclusive)
	theResults[3] == end matching index in reverse complement Seq
	theResults[4] == number of mismatches
*/
	var matLength=mat.length;
	var count=0;
	var theResults = new Array();
	var i;
	var	j;
	var	k;
	var mismatches;
  var maxInc;
	for(i=0; i<matLength; i++) {
		for(j=0; j<matLength; j++) {
			if(mat[i][j]==1)  { //potential start of an alignment
				mismatches=0;
				hasMatch=1;
				lastMatch=1;
				maxInc = matLength-(i<=j ? j : i);
				for (k=1; k<maxInc; k++) {
					hasMatch=mat[i+k][j+k];
					if(!hasMatch) break;
					if(hasMatch<=lastMatch) {
						if(mismatches>=maxMisMatch)
							break;
						mismatches++;
					}
					lastMatch=hasMatch;
				}
				if(k-mismatches>=minLen) {
					theResults[count]=new Array(5);
					theResults[count][0]=i;	//start index
					theResults[count][1]=j;	//start matching index in reverse complement seq
					theResults[count][2]=i+k-1; //end index of aligned bases (inclusive)
					theResults[count][3]=j+k-1; //end matching index in reverse complement Seq
					theResults[count][4]=mismatches;  //mismatch counts
					count++;
				}
			}
		}
	}
	return theResults;
}

function print_2d_string_array(array) {
    document.writeln("<table border>");

    var col_len = array.length;
    var row_len = array[0].length;

    for (var i = 0; i < row_len; i++) {
        document.writeln("<tr>");
        for (var j = 0; j < col_len; j++) {
           document.writeln("<td>" + array[j][i] + "</td>");
        }
        document.writeln("<tr>");
    }
    document.writeln("</table>");
}



var theSequence='GATCTCTGCAGAAGCGGAGC';
var theComplement = MakeComplement(theSequence, 1);
var theSeqArray=stringToArray(theSequence);
var theComplementArray=stringToArray(theComplement);

alert(theSequence);


var minAlignLen=5;
var maxMismatchNum=1;

var matrix= makeMatrix(theSequence.length);

fillMatchMatrix(theSeqArray, theComplementArray, matrix);

var mat=matrix;
print_2d_string_array(mat);

theAlignedArray = makeAlignedArray(matrix, minAlignLen, maxMismatchNum);
  
alert(theAlignedArray)
/*
alter(matrix)

theAlignedArray = makeAlignedArray(matrix, minAlignLen, maxMismatchNum);

*/
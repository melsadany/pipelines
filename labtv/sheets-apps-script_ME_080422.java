function createOneSlidePerRow() {


  let masterDeckID = "1zeQ7zIdpu4Ydcl5Ody_7uWxRGhbdrfl2LXc0ZjVTbxk";

  // Open the presentation and get the slides in it.
  let deck = SlidesApp.openById(masterDeckID);
  let slides = deck.getSlides();

  // The 2nd slide is the template that will be duplicated
  // once per row in the spreadsheet.
  let masterSlide = slides[1];
  for (var i = 2; i < slides.length; i++){
    slides[i].remove();
  };
  // let blankSlide = slides[2];


  // Load data from the spreadsheet.
  let dataRange = SpreadsheetApp.getActive().getDataRange();
  let sheetContents = dataRange.getValues();

  // Save the header in a variable called header
  let header = sheetContents.shift();

  // Create an array to save the data to be written back to the sheet.
  // I'll use this array to save links to the slides that are created.
  let updatedContents = [];

  // Reverse the order of rows because new slides will
  // be inserted at the top. Without this, the order of slides
  // will be the inverse of the ordering of rows in the sheet. 
  sheetContents.reverse();

  // For every row, create a new slide by duplicating the master slide
  // and replace the template variables with data from that row.
  sheetContents.forEach(function (row) {

    // Insert a new slide by duplicating the master slide.
    let slide = masterSlide.duplicate();

    // Populate data in the slide that was created
    slide.replaceAllText("{{title}}", row[0]);
    slide.replaceAllText("{{jabbrv}}", row[1]);
    slide.replaceAllText("{{month}}", row[2]);
    slide.replaceAllText("{{year}}", row[3]);
    slide.replaceAllText("{{abstract}}", row[4]);
    slide.replaceAllText("{{pmid}}", row[5]);
    

    updatedContents.push(row);
    
    // let imageSlide = blankSlide.duplicate();
    // imageSlide = deck.insertSlide(3)
    // imageSlide.insertImage(row[4]);
  });

  // Add the header back (remember it was removed using 
  // sheetContents.shift())
  updatedContents.push(header);

  // Reverse the array to preserve the original ordering of 
  // rows in the sheet.
  updatedContents.reverse();

  // Write the updated data back to the Google Sheets spreadsheet.
  dataRange.setValues(updatedContents);

  // Remove the master slide if you no longer need it.
  // masterSlide.remove();

}

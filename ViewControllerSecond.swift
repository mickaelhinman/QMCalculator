//
//  ViewControllerSecond.swift
//  calcGUI2
//
//  Created by (s) Mickael Hinman on 18/02/2016.
//  Copyright (c) 2016 (s) Mickael Hinman. All rights reserved.
//

import UIKit

class ViewControllerSecond: UIViewController {
    
    @IBOutlet weak var displaySecond: UILabel!
    
    var displayValue: String?
    
    var valuesToPlot: [[Float]]?
    
    @IBOutlet weak var graphPlot: Graph!
    
    @IBOutlet weak var lineThickness: UISlider!
    
    @IBAction func lineThicknessSlider(sender: AnyObject) {
        graphPlot.plotLineThickness = lineThickness.value
        graphPlot.setNeedsDisplay()
    }

    override func viewDidLoad() {
        super.viewDidLoad()
        
        displaySecond.text = displayValue
        print(valuesToPlot!)
        graphPlot.preparePointsToPlot(valuesToPlot!)
    }

    override func didReceiveMemoryWarning() {
        super.didReceiveMemoryWarning()
        // Dispose of any resources that can be recreated.
    }

}

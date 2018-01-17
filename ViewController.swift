//
//  ViewController.swift
//  calcGUI2
//
//  Created by (s) Mickael Hinman on 17/02/2016.
//  Copyright (c) 2016 (s) Mickael Hinman. All rights reserved.
//

import UIKit

class ViewController: UIViewController {
    
    var calculationArray: [EquationElement] = []
    
    @IBOutlet weak var basisStates: UITextField!
    
    @IBOutlet weak var display: UILabel!
    
    var displayText: [String] = []
    
    func addToDisplay(text: String){
        displayText.append(text)
        updateDisplay()
    }
    
    func updateDisplay(){
        var tempText = ""
        for(var index = 0; index < displayText.count; index += 1){
            tempText += displayText[index]
        }
        display.text = tempText
    }
    
    
    @IBAction func annihilationOperatorButton(sender: AnyObject) {
        addToDisplay("a ")
        calculationArray.append(annihilationOperator(basisStates.text.toInt()!))
    }
    
    @IBAction func creationOperatorButton(sender: AnyObject) {
        addToDisplay("aD ")
        calculationArray.append(creationOperator(basisStates.text.toInt()!))
    }
    
    @IBAction func numberOperatorButton(sender: AnyObject) {
        addToDisplay("n ")
        calculationArray.append(numberOperator(basisStates.text.toInt()!))
    }
    
    @IBAction func momentumOperatorpButton(sender: AnyObject) {
        addToDisplay("p ")
        calculationArray.append(momentumOperator(basisStates.text.toInt()!))
    }
    
    @IBAction func momentumOperatorPButton(sender: AnyObject) {
        addToDisplay("P ")
        calculationArray.append(momentumOperatorP(basisStates.text.toInt()!))
    }
    
    @IBAction func positionOperatorqButton(sender: AnyObject) {
        addToDisplay("q ")
        calculationArray.append(positionOperator(basisStates.text.toInt()!))
    }
    
    @IBAction func positionOperatorQButton(sender: AnyObject) {
        addToDisplay("Q ")
        calculationArray.append(positionOperatorQ(basisStates.text.toInt()!))
    }
    
    @IBAction func hButton(sender: AnyObject) {
        addToDisplay("h ")
        calculationArray.append(h)
    }
    
    @IBAction func hBarButton(sender: AnyObject) {
        addToDisplay("hBar ")
        calculationArray.append(hBar)
    }
    
    @IBAction func piButton(sender: AnyObject) {
        addToDisplay("pi ")
        calculationArray.append(pi)
    }
    
    @IBAction func iButton(sender: AnyObject) {
        addToDisplay("i ")
        calculationArray.append(i)
    }
    
    
    
    @IBAction func addButton(sender: AnyObject) {
        addToDisplay("+ ")
        calculationArray.append(EquationElement(.ADD))
    }
    
    @IBAction func subtractButton(sender: AnyObject) {
        addToDisplay("- ")
        calculationArray.append(EquationElement(.SUBTRACT))
    }
    
    @IBAction func divideButton(sender: AnyObject) {
        addToDisplay("/ ")
        calculationArray.append(EquationElement(.DIVIDE))
    }

    @IBAction func multiplyButton(sender: AnyObject) {
        addToDisplay("* ")
        calculationArray.append(EquationElement(.MULTIPLY))
    }
    
    @IBAction func transposeButton(sender: AnyObject) {
        addToDisplay("T ")
        calculationArray.append(EquationElement(.TRANSPOSE))
    }
    
    @IBAction func daggerButton(sender: AnyObject) {
        addToDisplay("D ")
        calculationArray.append(EquationElement(.DAGGER))
    }
    
    @IBAction func complexConjugateButton(sender: AnyObject) {
        addToDisplay("CC ")
        calculationArray.append(EquationElement(.COMPLEXCONJUGATE))
    }
    
    @IBAction func negateButton(sender: AnyObject) {
        addToDisplay("- ")
        calculationArray.append(EquationElement(.NEGATE))
    }
    
    @IBAction func squaredButton(sender: AnyObject) {
        addToDisplay("^2 ")
        calculationArray.append(EquationElement(.SQUARED))
    }
    
    @IBAction func cubedButton(sender: AnyObject) {
        addToDisplay("^3 ")
        calculationArray.append(EquationElement(.CUBED))
    }
    
    @IBAction func powerButton(sender: AnyObject) {
        addToDisplay("^ ")
        calculationArray.append(EquationElement(.POWER))
    }
    
    @IBAction func squareRootButton(sender: AnyObject) {
        addToDisplay("root ")
        calculationArray.append(EquationElement(.SQUAREROOT))
    }
    
    @IBAction func factorialButton(sender: AnyObject) {
        addToDisplay("! ")
        calculationArray.append(EquationElement(.FACTORIAL))
    }
    
    @IBAction func leftBracketButton(sender: AnyObject) {
        addToDisplay("( ")
        calculationArray.append(EquationElement(.LEFT))
    }
    
    @IBAction func rightBracketButton(sender: AnyObject) {
        addToDisplay(") ")
        calculationArray.append(EquationElement(.RIGHT))
    }
    
    @IBAction func sin(sender: AnyObject) {
        addToDisplay("sin ")
        calculationArray.append(EquationElement(.SIN))
    }
    
    @IBAction func cos(sender: AnyObject) {
        addToDisplay("cos ")
        calculationArray.append(EquationElement(.COS))
    }
    
    @IBAction func exp(sender: AnyObject) {
        addToDisplay("exp ")
        calculationArray.append(EquationElement(.EXP))
    }
    
    @IBAction func identity(sender: AnyObject) {
        addToDisplay("Id ")
        calculationArray.append(identityMatrix(basisStates.text.toInt()!))
    }
    
    var userConstantNumber = 0
    
    func getUserConstantNumber() -> Int{
        userConstantNumber += 1
        return userConstantNumber
    }
    
    @IBOutlet weak var userConstant: UITextField!
    
    @IBAction func matrixButton(sender: AnyObject) {
        let values = userConstant.text.componentsSeparatedByString(",")
        var valuesFloatArray: [Float] = []
        for(var index = 0; index < values.count; index += 1){
            valuesFloatArray.append(NSNumberFormatter().numberFromString(values[index])!.floatValue)
        }
        calculationArray.append(EquationElement(.MATRIX,valuesFloatArray))
        addToDisplay("M" + String(stringInterpolationSegment: getUserConstantNumber()) + " ")
        userConstant.text = ""
    }
    
    @IBAction func matrixComplexButton(sender: AnyObject) {
        let values = userConstant.text.componentsSeparatedByString(",")
        var valuesFloatArray: [Float] = []
        for(var index = 0; index < values.count; index += 1){
            valuesFloatArray.append(NSNumberFormatter().numberFromString(values[index])!.floatValue)
        }
        calculationArray.append(EquationElement(.MATRIXCOMPLEX,valuesFloatArray))
        addToDisplay("MZ" + String(stringInterpolationSegment: getUserConstantNumber()) + " ")
        userConstant.text = ""
    }
    
    @IBAction func vectorButton(sender: AnyObject) {
        let values = userConstant.text.componentsSeparatedByString(",")
        var valuesFloatArray: [Float] = []
        for(var index = 0; index < values.count; index += 1){
            valuesFloatArray.append(NSNumberFormatter().numberFromString(values[index])!.floatValue)
        }
        calculationArray.append(EquationElement(.VECTOR,valuesFloatArray))
        addToDisplay("V" + String(stringInterpolationSegment: getUserConstantNumber()) + " ")
        userConstant.text = ""
    }
    
    @IBAction func vectorComplexButton(sender: AnyObject) {
        let values = userConstant.text.componentsSeparatedByString(",")
        var valuesFloatArray: [Float] = []
        for(var index = 0; index < values.count; index += 1){
            valuesFloatArray.append(NSNumberFormatter().numberFromString(values[index])!.floatValue)
        }
        calculationArray.append(EquationElement(.VECTORCOMPLEX,valuesFloatArray))
        addToDisplay("VZ" + String(stringInterpolationSegment: getUserConstantNumber()) + " ")
        userConstant.text = ""
    }
    
    @IBAction func constantButton(sender: AnyObject) {
        let value = userConstant.text
        let valueFloat: Float
        valueFloat = NSNumberFormatter().numberFromString(value)!.floatValue
        
        calculationArray.append(EquationElement(valueFloat))
        addToDisplay("C" + String(stringInterpolationSegment: getUserConstantNumber()) + " ")
        userConstant.text = ""
    }
    
    @IBAction func constantComplexButton(sender: AnyObject) {
        let values = userConstant.text.componentsSeparatedByString(",")
        var valuesFloatArray: [Float] = []
        for(var index = 0; index < values.count; index += 1){
            valuesFloatArray.append(NSNumberFormatter().numberFromString(values[index])!.floatValue)
        }
        calculationArray.append(EquationElement(.CONSTANTCOMPLEX,valuesFloatArray))
        addToDisplay("CZ" + String(stringInterpolationSegment: getUserConstantNumber()) + " ")
        userConstant.text = ""
    }
    
    
    
    
    @IBOutlet weak var startValue: UITextField!
    
    @IBOutlet weak var incrementValue: UITextField!
    
    @IBOutlet weak var numberOfSteps: UITextField!
    
    var variableIndex = -1
    var plotStartValue: Float = -1.0
    var plotIncrementValue: Float = -1.0
    var plotNumberOfSteps: Float = -1.0
        
    @IBAction func addVariableButton(sender: AnyObject) {
        
        plotStartValue = NSNumberFormatter().numberFromString(startValue.text)!.floatValue
        plotIncrementValue = NSNumberFormatter().numberFromString(incrementValue.text)!.floatValue
        plotNumberOfSteps = NSNumberFormatter().numberFromString(numberOfSteps.text)!.floatValue
        
        variableIndex = calculationArray.count
        addToDisplay("var ")
    }
    
    
    @IBAction func backspaceButton(sender: AnyObject) {
        calculationArray.removeLast()
        displayText.removeLast()
        updateDisplay()
    }
    
    @IBAction func calculateButton(sender: AnyObject) {
        let result = implimentElementArraySimplification(calculationArray)
        var tempResult = ""
        
        if(result.type == .CONSTANT){
            tempResult += "Constant -> "
            tempResult += String(stringInterpolationSegment: result.constant)
        }
        else if(result.type == .CONSTANTCOMPLEX){
            tempResult += "Constant Complex -> "
            for(var index = 0; index < result.constantComplex.count; index += 1){
                tempResult += String(stringInterpolationSegment: result.constantComplex[index]) + ","
            }
        }
        else if(result.type == .VECTOR){
            tempResult += "Vector -> "
            for(var index = 0; index < result.vector.count; index += 1){
                tempResult += String(stringInterpolationSegment: result.vector[index]) + ","
            }
        }
        else if(result.type == .VECTORCOMPLEX){
            tempResult += "Vector Complex -> "
            for(var index = 0; index < result.vectorComplex.count; index += 1){
                tempResult += String(stringInterpolationSegment: result.vectorComplex[index]) + ","
            }
        }
        else if(result.type == .MATRIX){
            tempResult += "Matrix -> "
            for(var index = 0; index < result.matrix.count; index += 1){
                tempResult += String(stringInterpolationSegment: result.matrix[index]) + ","
            }
        }
        else if(result.type == .MATRIXCOMPLEX){
            tempResult += "Matrix Complex -> "
            for(var index = 0; index < result.matrixComplex.count; index += 1){
                tempResult += String(stringInterpolationSegment: result.matrixComplex[index]) + ","
            }
        }
        display.text = tempResult
    }
    
    
    
    override func prepareForSegue(segue: UIStoryboardSegue, sender: AnyObject?) {
        let destination = segue.destinationViewController as! ViewControllerSecond
        
        updateDisplay()
        
        destination.displayValue = display.text
        destination.valuesToPlot = plotableEigenValues(calculationArray, plotStartValue, plotIncrementValue, plotNumberOfSteps, variableIndex)
    }
    
    override func viewDidLoad() {
        super.viewDidLoad()
        // Do any additional setup after loading the view, typically from a nib.
    }

    override func didReceiveMemoryWarning() {
        super.didReceiveMemoryWarning()
        // Dispose of any resources that can be recreated.
    }


}


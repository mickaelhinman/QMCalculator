//
//  Graph.swift
//  calcGUI2
//
//  Created by (s) Mickael Hinman on 18/02/2016.
//  Copyright (c) 2016 (s) Mickael Hinman. All rights reserved.
//

import UIKit

class Graph: UIView {
    
    //let plotWidth = UIScreen.mainScreen().bounds.size.width
    let plotWidth: Float = 730
    let plotHeight: Float = 450
    
    var xAxisOrigin: Float = 0
    
    var finalPoints: [[Float]] = []
    
    func findMostNegativePoint(points: [Float]) -> Float{
        var mostNegativePoint: Float = 0
        for(var index = 0; index < points.count; index += 1){
            if(points[index] < mostNegativePoint){
                mostNegativePoint = points[index]
            }
        }
        return -mostNegativePoint
    }
    
    func findLeastPositivePoint(points: [Float]) -> Float{
        var leastPositivePoint: Float = points[0]
        for(var index = 0; index < points.count; index += 1){
            if(points[index] < leastPositivePoint){
                leastPositivePoint = points[index]
            }
        }
        return leastPositivePoint
    }
    
    func findMostPositivePoint(points: [Float]) -> Float{
        var mostPositivePoint: Float = 1
        for(var index = 0; index < points.count; index += 1){
            if(points[index] > mostPositivePoint){
                mostPositivePoint = points[index]
            }
        }
        return mostPositivePoint
    }
    
    func preparePointsToPlot(var valuesToPlotArray: [[Float]]){
        var valuesToPlot: [Float] = []
        for(var index = 0; index < valuesToPlotArray.count; index += 1){
            for(var place = 0; place < valuesToPlotArray[0].count; place += 1){
                valuesToPlot.append(valuesToPlotArray[index][place])
            }
        }
        
        var addOnValue = findMostNegativePoint(valuesToPlot)
        
        if(addOnValue == 0){
            addOnValue = -1 * findLeastPositivePoint(valuesToPlot)
        }
        
        for(var index = 0; index < valuesToPlot.count; index += 1){
            valuesToPlot[index] += addOnValue
        }
        
        var scaleValue: Float
        scaleValue = plotHeight/findMostPositivePoint(valuesToPlot)
        
        
        for(var index = 0; index < valuesToPlot.count; index += 1){
            valuesToPlot[index] *= scaleValue
        }
        for(var index = 0; index < valuesToPlot.count; index += 1){
            valuesToPlot[index] = plotHeight - valuesToPlot[index]
        }
        
        
        var finalPointsToPlot: [[Float]] = []
        var binNumber = 0
        for(var index = 0; index < valuesToPlot.count; index += 1){
            if(index == 0){
                for(var place = 0; place < valuesToPlotArray.count; place += 1){
                    finalPointsToPlot.append([])
                }
            }
            if(index % valuesToPlotArray[0].count == 0 && index != 0){
                binNumber += 1
            }
            
            finalPointsToPlot[binNumber].append(valuesToPlot[index])
        }
        
        finalPoints = finalPointsToPlot
        print(finalPoints)
    }
    
    
    
    
    var plotLineThickness: Float = 2.0
    
    
    override func drawRect(rect: CGRect) {
        let context = UIGraphicsGetCurrentContext()
        CGContextSetLineWidth(context, CGFloat(plotLineThickness))
        
        let colorSpace = CGColorSpaceCreateDeviceRGB()
        //let components: [CGFloat] = [0.0, 0.0, 0.0, 1.0]
        //let color = CGColorCreate(colorSpace, components)
        //CGContextSetStrokeColorWithColor(context, color)
        
        
        //CGContextMoveToPoint(context, 30, 0)
        //CGContextAddLineToPoint(context, 30, 470)
        //CGContextMoveToPoint(context, 0, CGFloat(xAxisOrigin))
        //CGContextAddLineToPoint(context, screenWidth, CGFloat(xAxisOrigin))
        
        let componentsPoints: [CGFloat] = [0.0, 0.0, 0.0, 1.0]
        let colorPoints = CGColorCreate(colorSpace, componentsPoints)
        CGContextSetStrokeColorWithColor(context, colorPoints)
        
        let numberOfPoints = finalPoints[0].count
        let interval = Float(plotWidth) / Float(numberOfPoints - 1)
        
        for(var plotLine = 0; plotLine < finalPoints.count; plotLine += 1){
            for(var index = 0; index < numberOfPoints; index += 1){
                if(index == 0){
                    CGContextMoveToPoint(context, CGFloat(0), CGFloat(finalPoints[plotLine][index]))
                }
                else{
                    CGContextAddLineToPoint(context, CGFloat((Float(index) * interval) + 0), CGFloat(finalPoints[plotLine][index]))
                }
            }
        }
        
        
        CGContextStrokePath(context)
    }
    

}
